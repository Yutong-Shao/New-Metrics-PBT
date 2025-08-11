#!/usr/bin/env python3

import os
import sys
import math
import time
import glob
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
from scipy.stats import cramervonmises_2samp

printed_model_adequacy_header = False

# IO Utilities

def read_fasta(filepath):
    sequences = []
    with open(filepath, "r") as file:
        for line in file:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                sequences.append({"header": line, "sequence": ""})
            else:
                sequences[-1]["sequence"] += line
    return sequences

def write_fasta(filepath, records):
    with open(filepath, "w") as file:
        for rec in records:
            file.write(f"{rec['header']}\n{rec['sequence']}\n")

# Gap Transfer Logic

## Extract gap positions from the original alignment (by sequence index)
## We treat characters -, X, ? as gaps/missing

def get_gap_positions(alignment):
    gap_positions = []
    for entry in alignment:
        pos = [i for i, aa in enumerate(entry["sequence"]) if aa in {"-", "X", "?"}]
        gap_positions.append(pos)
    return gap_positions

def transfer_gaps(target_alignment, gap_positions):
    out = []
    for idx, entry in enumerate(target_alignment):
        seq = list(entry["sequence"])
        for p in gap_positions[idx]:
            if p < len(seq):
                seq[p] = "-"
        out.append({"header": entry["header"], "sequence": "".join(seq)})
    return out

def read_fasta_dict(filepath):
    with open(filepath, "r") as f:
        return {title: seq for title, seq in SimpleFastaParser(f)}

def stage1_add_gaps(sim_root_dir, original_alignment_path):
    print(f"\n[STEP 1] Loading original alignment: {original_alignment_path}")
    orig = read_fasta(original_alignment_path)
    headers = [r["header"] for r in orig]
    gap_pos = get_gap_positions(orig)

    for model in sorted(os.listdir(sim_root_dir)):
        mp = os.path.join(sim_root_dir, model)
        if not os.path.isdir(mp):
            continue
        print(f"[MODEL] Adding gaps in: {model}")
        for fn in sorted(os.listdir(mp)):
            if not (fn.endswith(".fa") or fn.endswith(".fasta")):
                continue
            path = os.path.join(mp, fn)
            sim = read_fasta(path)
            if len(sim) != len(orig):
                continue
            if any(len(sim[i]["sequence"]) != len(orig[i]["sequence"]) for i in range(len(orig))):
                continue
            for i in range(len(sim)):
                sim[i]["header"] = headers[i]
            gapped = transfer_gaps(sim, gap_pos)
            outfn = f"{fn}_gap"
            write_fasta(os.path.join(mp, outfn), gapped)

# ---------- Step 2: Mean entropy & diversity in parallel ----------

def calc_shannon_entropy(column):
    counts = Counter(aa for aa in column if aa not in {"-", "X", "?"})
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = -sum((count / total) * math.log2(count / total) for count in counts.values())
    return entropy

def calc_site_diversity(column):
    aas = [aa for aa in column if aa not in {"-", "X", "?"}]
    return len(set(aas))

def calc_mean_entropy(sequences):
    alignment_length = len(next(iter(sequences.values())))
    total_entropy = sum(
        calc_shannon_entropy([seq[i] for seq in sequences.values()])
        for i in range(alignment_length)
    )
    return total_entropy / alignment_length

def calc_mean_diversity(sequences):
    alignment_length = len(next(iter(sequences.values())))
    total_div = sum(
        calc_site_diversity([seq[i] for seq in sequences.values()])
        for i in range(alignment_length)
    )
    return total_div / alignment_length

def compute_zscore(original, simulated_list):
    avg = sum(simulated_list) / len(simulated_list)
    var = sum((x - avg) ** 2 for x in simulated_list) / (len(simulated_list) - 1)
    sd = math.sqrt(var)
    z = -(original - avg) / sd if sd != 0 else float("nan")
    return avg, sd, z

def process_model_folder(model_path, original_alignment):
    entropy_vals = []
    diversity_vals = []
    entropy_files = []
    diversity_files = []

    for file in os.listdir(model_path):
        if file.endswith(".fa_gap"):
            full_path = os.path.join(model_path, file)
            sequences = read_fasta_dict(full_path)
            entropy = calc_mean_entropy(sequences)
            diversity = calc_mean_diversity(sequences)
            entropy_vals.append(entropy)
            entropy_files.append(file)
            diversity_vals.append(diversity)
            diversity_files.append(file)

    if entropy_vals:
        avg_e, sd_e, z_e = compute_zscore(calc_mean_entropy(original_alignment), entropy_vals)
        with open(os.path.join(model_path, "entropy.pbr_gaps"), "w") as f:
            f.write("Test of model adequacy / across sites entropy heterogeneity\n")
            f.write(f"Entropy (original data): {calc_mean_entropy(original_alignment)}\n")
            f.write(f"Average entropy (simulated data): {avg_e}\n")
            f.write(f"SD simulated data: {sd_e}\n")
            f.write(f"Z-score: {z_e}\n")

        with open(os.path.join(model_path, "entropy_scores_bootstrapped_data.txt_gaps"), "w") as f:
            f.write("file\tentropy\n")
            for file, val in zip(entropy_files, entropy_vals):
                f.write(f"{file}\t{val}\n")

    if diversity_vals:
        avg_d, sd_d, z_d = compute_zscore(calc_mean_diversity(original_alignment), diversity_vals)
        with open(os.path.join(model_path, "diversity.pbr_gaps"), "w") as f:
            f.write("Test of model adequacy / across sites compositional heterogeneity\n")
            f.write(f"Diversity original data: {calc_mean_diversity(original_alignment)}\n")
            f.write(f"Average diversity simulated data: {avg_d}\n")
            f.write(f"SD simulated data: {sd_d}\n")
            f.write(f"Z-score: {z_d}\n")

        with open(os.path.join(model_path, "diversity_scores_bootstrapped_data.txt_gaps"), "w") as f:
            f.write("file\tdiversity\n")
            for file, val in zip(diversity_files, diversity_vals):
                f.write(f"{file}\t{val}\n")

def stage2_parallel_entropy_diversity(sim_root, orig_path, num_cpus):
    original_alignment = read_fasta_dict(orig_path)
    model_dirs = [os.path.join(sim_root, d) for d in sorted(os.listdir(sim_root)) if os.path.isdir(os.path.join(sim_root, d))]
    print(f"\n[STEP 2] Computing entropy & diversity using {num_cpus} CPUs...")
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        for model_dir in model_dirs:
            executor.submit(process_model_folder, model_dir, original_alignment)

def summarize_mean_values(sim_root_dir):
    diversity_summary = []
    entropy_summary = []

    for model in sorted(os.listdir(sim_root_dir)):
        model_path = os.path.join(sim_root_dir, model)
        if not os.path.isdir(model_path):
            continue

        # Diversity
        div_path = os.path.join(model_path, "diversity.pbr_gaps")
        if os.path.isfile(div_path):
            with open(div_path, "r") as f:
                lines = f.readlines()
                try:
                    div_orig = float(lines[1].split(":")[1].strip())
                    div_avg = float(lines[2].split(":")[1].strip())
                    div_sd = float(lines[3].split(":")[1].strip())
                    div_z = float(lines[4].split(":")[1].strip())
                    diversity_summary.append([model, div_orig, div_avg, div_sd, div_z])
                except:
                    print(f"[WARNING] Failed to parse {div_path}")

        # Entropy
        ent_path = os.path.join(model_path, "entropy.pbr_gaps")
        if os.path.isfile(ent_path):
            with open(ent_path, "r") as f:
                lines = f.readlines()
                try:
                    ent_orig = float(lines[1].split(":")[1].strip())
                    ent_avg = float(lines[2].split(":")[1].strip())
                    ent_sd = float(lines[3].split(":")[1].strip())
                    ent_z = float(lines[4].split(":")[1].strip())
                    entropy_summary.append([model, ent_orig, ent_avg, ent_sd, ent_z])
                except:
                    print(f"[WARNING] Failed to parse {ent_path}")

    # Write diversity summary
    with open(os.path.join(sim_root_dir, "mean_diversity_results.txt"), "w") as f:
        f.write("Model\tDiversity original data\tAverage diversity simulated data\tSD simulated data\tZ-score\n")
        for row in diversity_summary:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\n")

    # Write entropy summary
    with open(os.path.join(sim_root_dir, "mean_entropy_results.txt"), "w") as f:
        f.write("Model\tEntropy original data\tAverage entropy simulated data\tSD simulated data\tZ-score\n")
        for row in entropy_summary:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\n")


# ---------- Step 3a: Original per-site entropy & diversity ----------

def compute_original_sitewise(sim_root, orig_path):
    print("\n[STEP 3a] Computing per-site entropy & diversity for original alignment")
    with open(orig_path, "r") as f:
        seqs = {h: s for h, s in SimpleFastaParser(f)}
    L = len(next(iter(seqs.values())))

    entropy_vals = [calc_shannon_entropy([s[i] for s in seqs.values()]) for i in range(L)]
    with open(os.path.join(sim_root, "original.sitewise_entropy.txt"), "w") as f:
        f.write("site\tentropy\n")
        for i, v in enumerate(entropy_vals, 1):
            f.write(f"{i}\t{v:.6f}\n")

    diversity_vals = [calc_site_diversity([s[i] for s in seqs.values()]) for i in range(L)]
    with open(os.path.join(sim_root, "original.sitewise_diversity.txt"), "w") as f:
        f.write("site\tdiversity\n")
        for i, v in enumerate(diversity_vals, 1):
            f.write(f"{i}\t{v}\n")

# ---------- Step 3b: Per-site metrics in parallel ----------

def site_diversity(column):
    aas = [aa for aa in column if aa not in {"-", "X", "?"}]
    return len(set(aas))

def site_entropy(column):
    cnt = Counter(aa for aa in column if aa not in {"-", "X", "?"})
    tot = sum(cnt.values())
    if tot == 0:
        return 0.0
    h = 0.0
    for v in cnt.values():
        p = v / tot
        h -= p * math.log2(p)
    return h

def compute_sitewise_diversity(fasta_path):
    with open(fasta_path, "r") as f:
        seqs = {h: s for h, s in SimpleFastaParser(f)}
    L = len(next(iter(seqs.values())))
    divs = [site_diversity([s[i] for s in seqs.values()]) for i in range(L)]
    with open(fasta_path + ".sitewise_diversity.txt", "w") as out:
        out.write("site\tdiversity\n")
        for i, v in enumerate(divs, 1):
            out.write(f"{i}\t{v}\n")

def compute_sitewise_entropy(fasta_path):
    with open(fasta_path, "r") as f:
        seqs = {h: s for h, s in SimpleFastaParser(f)}
    L = len(next(iter(seqs.values())))
    ents = [site_entropy([s[i] for s in seqs.values()]) for i in range(L)]
    with open(fasta_path + ".sitewise_entropy.txt", "w") as out:
        out.write("site\tentropy\n")
        for i, v in enumerate(ents, 1):
            out.write(f"{i}\t{v:.6f}\n")

def process_gap_file(path):
    compute_sitewise_diversity(path)
    compute_sitewise_entropy(path)

def stage3_compute_metrics_parallel(sim_root_dir, num_cpus):
    gap_files = []
    for model in sorted(os.listdir(sim_root_dir)):
        mp = os.path.join(sim_root_dir, model)
        if not os.path.isdir(mp):
            continue
        for fn in os.listdir(mp):
            if fn.endswith(".fa_gap"):
                gap_files.append(os.path.join(mp, fn))

    print(f"\n[STEP 3b] Computing per-site metrics on {len(gap_files)} files using {num_cpus} CPUs...")
    with ProcessPoolExecutor(max_workers=num_cpus) as ex:
        ex.map(process_gap_file, gap_files)


# ---------- Step 4: CvM Test ----------

def cvm_test_per_model(sim_root_dir):
    print("\n[STEP 4] Running CvM test comparing model sitewise metrics to original")
    # Load original
    orig_div = pd.read_csv(os.path.join(sim_root_dir, "original.sitewise_diversity.txt"), sep='\t')["diversity"].values
    orig_ent = pd.read_csv(os.path.join(sim_root_dir, "original.sitewise_entropy.txt"), sep='\t')["entropy"].values

    entropy_results = []
    diversity_results = []

    for model in sorted(os.listdir(sim_root_dir)):
        model_path = os.path.join(sim_root_dir, model)
        if not os.path.isdir(model_path):
            continue
        for file in os.listdir(model_path):
            if file.endswith(".fa_gap.sitewise_entropy.txt"):
                df = pd.read_csv(os.path.join(model_path, file), sep='\t')
                seq_vals = df["entropy"].values
                result = cramervonmises_2samp(seq_vals, orig_ent)
                entropy_results.append((model, file, result.statistic, result.pvalue))
            elif file.endswith(".fa_gap.sitewise_diversity.txt"):
                df = pd.read_csv(os.path.join(model_path, file), sep='\t')
                seq_vals = df["diversity"].values
                result = cramervonmises_2samp(seq_vals, orig_div)
                diversity_results.append((model, file, result.statistic, result.pvalue))

    with open(os.path.join(sim_root_dir, "cvm_entropy_results.txt"), "w") as f:
        f.write("model\tfile\tW2_statistic\tp_value\n")
        for m, f_, s, p in entropy_results:
            f.write(f"{m}\t{f_}\t{s:.6e}\t{p:.6e}\n")

    with open(os.path.join(sim_root_dir, "cvm_diversity_results.txt"), "w") as f:
        f.write("model\tfile\tW2_statistic\tp_value\n")
        for m, f_, s, p in diversity_results:
            f.write(f"{m}\t{f_}\t{s:.6e}\t{p:.6e}\n")


# ---------- Step 5: Find Best-fit Model and Print Sorted List ----------

def find_best_fit_model_by_mean_difference(file_path, metric_name):
    global printed_model_adequacy_header

    df = pd.read_csv(file_path, sep="\t")

    if metric_name.lower() == "diversity":
        orig_col = "Diversity original data"
        avg_col = "Average diversity simulated data"
    elif metric_name.lower() == "entropy":
        orig_col = "Entropy original data"
        avg_col = "Average entropy simulated data"
    else:
        print(f"[ERROR] Unknown metric: {metric_name}")
        return

    df["abs_diff"] = (df[orig_col] - df[avg_col]).abs()
    df_sorted = df.sort_values(by="abs_diff")

    best_model = df_sorted.iloc[0]["Model"]

    if not printed_model_adequacy_header:
        print()
        print()
        print("Model Adequacy Evaluation")
        print("----------------------")
        printed_model_adequacy_header = True  # set flag so it won't print again

    print()
    print()
    print(f"Best-fit model according to Mean Difference ({metric_name}): {best_model}")
    print()
    print(f"List of models sorted by |original - simulated| for {metric_name}:")
    print()
    header = f"{'Model':<35} {orig_col:<25} {avg_col:<30} {'SD simulated':<18} {'Z-score':<10} {'|Diff|':<10}"
    print(header)
    print("-" * len(header))
    print()

    for _, row in df_sorted.iterrows():
        print(f"{row['Model']:<35} {row[orig_col]:<25.6f} {row[avg_col]:<30.6f} {row['SD simulated data']:<18.6f} {row['Z-score']:<10.3f} {row['abs_diff']:<10.6f}")

def find_best_fit_model_by_cvm(file_path, metric_name):
    global printed_model_adequacy_header

    df = pd.read_csv(file_path, sep="\t")
    model_means = df.groupby("model")["W2_statistic"].mean()
    sorted_models = model_means.sort_values()

    best_model = sorted_models.index[0]

    if not printed_model_adequacy_header:
        print()
        print()
        print("Model Adequacy Evaluation")
        print("----------------------")
        printed_model_adequacy_header = True  # set flag so it won't print again

    print()
    print()
    print(f"Best-fit model according to Cvm Test ({metric_name}): {best_model}")
    print()
    print(f"List of models sorted by Cvm Test ({metric_name}) W2_statistic:")
    print()
    print(f"{'Model':<35} {'W2_statistic':>15}")
    for model, val in sorted_models.items():
        print(f"{model:<35} {val:>15.6f}")
    print()

# Main

def main():
    if len(sys.argv) < 3:
        print("Usage: python PBT.py <SimRootFolder> <OriginalAlignment> [NumCPUs]")
        sys.exit(1)

    sim_root = sys.argv[1]
    orig_path = sys.argv[2]
    cpus = int(sys.argv[3]) if len(sys.argv) > 3 else 1

    if not os.path.isdir(sim_root):
        print(f"[ERROR] Simulation folder not found: {sim_root}")
        sys.exit(1)
    if not os.path.isfile(orig_path):
        print(f"[ERROR] Original alignment not found: {orig_path}")
        sys.exit(1)

    # Start time
    start_time = time.time()

    # Logging setup
    class Logger(object):
        def __init__(self, logfile_path):
            self.terminal = sys.stdout
            self.log = open(logfile_path, "w", encoding="utf-8")
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
        def flush(self):
            self.terminal.flush()
            self.log.flush()

    log_path = os.path.join(sim_root, "PBT.log")
    sys.stdout = Logger(log_path)

    # Log header
    print("=" * 60)
    print(f"DATE AND TIME: {time.strftime('%a %b %d %H:%M:%S %Y', time.localtime(start_time))}")
    print("Starting parametric bootstrap test for model adequacy")
    print(f"Simulated folder path: {sim_root}")
    print(f"Original alignment path: {orig_path}")
    print("=" * 60 + "\n")

    # Actual steps
    stage1_add_gaps(sim_root, orig_path)
    stage2_parallel_entropy_diversity(sim_root, orig_path, cpus)
    summarize_mean_values(sim_root)
    compute_original_sitewise(sim_root, orig_path)
    stage3_compute_metrics_parallel(sim_root, cpus)
    cvm_test_per_model(sim_root)

    # Step 5: Print mean-based best-fit model and sorted list
    mean_entropy_file = os.path.join(sim_root, "mean_entropy_results.txt")
    mean_diversity_file = os.path.join(sim_root, "mean_diversity_results.txt")
    if os.path.isfile(mean_entropy_file):
        find_best_fit_model_by_mean_difference(mean_entropy_file, "entropy")
    if os.path.isfile(mean_diversity_file):
        find_best_fit_model_by_mean_difference(mean_diversity_file, "diversity")

    # Step 6: Print best-fit model and sorted list
    entropy_file = os.path.join(sim_root, "cvm_entropy_results.txt")
    diversity_file = os.path.join(sim_root, "cvm_diversity_results.txt")
    if os.path.isfile(entropy_file):
        find_best_fit_model_by_cvm(entropy_file, "entropy")
    if os.path.isfile(diversity_file):
        find_best_fit_model_by_cvm(diversity_file, "diversity")

    # Final log time
    end_time = time.time()
    elapsed_time = end_time - start_time
    timestamp = time.strftime("%a %b %d %H:%M:%S %Y", time.localtime(end_time))

    print("\nTIME STAMP")
    print("----------")
    print(f"Date and time: {timestamp}")
    print(f"Total wall-clock time used: {elapsed_time:.5f} seconds ({int(elapsed_time//3600)}h:{int((elapsed_time%3600)//60)}m:{int(elapsed_time%60)}s)")

if __name__ == "__main__":
    main()

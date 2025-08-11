import os
import sys
import csv
import re
from collections import defaultdict

def clean_rate_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    first_nine = lines[:9]
    should_clean = True
    for line in first_nine:
        line = line.strip()
        if not (line.startswith("#") or line == "Site\tRate\tCat\tC_Rate"):
            should_clean = False
            break

    if should_clean:
        new_lines = lines[9:]
        with open(file_path, 'w') as f:
            f.writelines(new_lines)
        print(f"[CLEANED] {file_path}")
    else:
        print(f"[SKIPPED] {file_path}")

def convert_sitefreq(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    formatted_lines = [" ".join(line.strip().split()) + "\n" for line in lines]

    with open(file_path, 'w', encoding='utf-8') as file:
        file.writelines(formatted_lines)

    print(f"[REFORMATTED] {file_path}")

def parse_iqtree_model_info(iqtree_file):
    with open(iqtree_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    model_prefix = "UnknownModel"
    model_suffix = "UNKNOWN"
    site_rates = "{}"

    for line in lines:
        if line.startswith("Model of substitution:"):
            model_line = line.replace("Model of substitution:", "").strip()
            parts = model_line.split("+")
            if len(parts) >= 2:
                model_prefix = parts[0]
                model_suffix = parts[-1]

        elif line.startswith("Site proportion and rates:"):
            values = re.findall(r"[\d.]+", line)
            site_rates = "{" + ",".join(values) + "}"

    return model_prefix, model_suffix, site_rates

def create_nexus(sitefreq_file, rate_file, iqtree_file, output_file):
    my_dict = defaultdict(list)

    with open(sitefreq_file, 'r') as frequencies:
        reader = csv.reader(frequencies, delimiter=' ')
        for row in reader:
            my_dict[row[0]].append(row[1:])

    with open(rate_file, 'r') as rates:
        reader = csv.reader(rates, delimiter='\t')
        for row in reader:
            site = row[0]
            rate = float(row[1])
            gamma = float(row[3])
            my_dict[site].append(rate)
            my_dict[site].append(gamma)

    model_prefix, model_suffix, site_rates = parse_iqtree_model_info(iqtree_file)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('#nexus\nbegin sets;\n')

        for k in my_dict:
            f.write(f"    charset site_{k} = {k};\n")

        f.write('    charpartition mine = \n')

        items = list(my_dict.items())
        for idx, (k, v) in enumerate(items):
            freq = '/'.join(v[0])
            tree_len = str(v[1])
            gamma = str(v[2])

            if model_suffix.startswith("R"):
                model_str = f"{model_prefix}+F{{{freq}}}+{model_suffix}{site_rates}"
            elif model_suffix.startswith("G"):
                model_str = f"{model_prefix}+F{{{freq}}}+G{{{gamma}}}"
            else:
                model_str = f"{model_prefix}+F{{{freq}}}+UNKNOWN"

            end = ";\n" if idx == len(items) - 1 else ",\n"
            f.write(f"        {model_str}:site_{k}{{{tree_len}}}{end}")

        f.write('end;\n')
    print(f"[CREATED] {output_file}")

def convert_nex_format(file_path):
    with open(file_path, "r", encoding="utf-8") as file:
        content = file.read()

    updated_content = re.sub(r'(F\{[^}]+)', lambda m: m.group(1).replace('/', ','), content)

    with open(file_path, "w", encoding="utf-8") as file:
        file.write(updated_content)

    print(f"[FORMATTED] {file_path}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python combine_sitefq_tate.py <input_folder>")
        sys.exit(1)

    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"[ERROR] Provided path is not a valid directory: {input_dir}")
        sys.exit(1)

    files = os.listdir(input_dir)
    sitefreq_files = [f for f in files if f.endswith("step1.sitefreq")]
    rate_files = [f for f in files if f.endswith("step2.rate")]
    iqtree_files = [f for f in files if f.endswith("step2.iqtree")]

    sitefreq_map = {f.replace("-step1.sitefreq", ""): os.path.join(input_dir, f) for f in sitefreq_files}
    rate_map = {f.replace("-step2.rate", ""): os.path.join(input_dir, f) for f in rate_files}
    iqtree_map = {f.replace("-step2.iqtree", ""): os.path.join(input_dir, f) for f in iqtree_files}

    print("\nStep 1: Cleaning .rate files")
    for rate_file in rate_map.values():
        clean_rate_file(rate_file)

    print("\nStep 2: Formatting .sitefreq files")
    for sitefreq_file in sitefreq_map.values():
        convert_sitefreq(sitefreq_file)

    print("\nStep 3: Creating .nex files")
    for model in sitefreq_map:
        if model in rate_map and model in iqtree_map:
            sitefreq_file = sitefreq_map[model]
            rate_file = rate_map[model]
            iqtree_file = iqtree_map[model]
            output_file = os.path.join(input_dir, f"{model}.nex")

            print(f"\nGenerating NEXUS file for model: {model}")
            create_nexus(sitefreq_file, rate_file, iqtree_file, output_file)
            convert_nex_format(output_file)
        else:
            print(f"\n[WARNING] Missing file(s) for model: {model}")

if __name__ == "__main__":
    main()

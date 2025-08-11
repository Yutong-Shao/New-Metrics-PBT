import os
import re
import sys

def extract_fmix_from_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    fmix_entries = []
    capture = False

    for line in lines:
        if re.match(r"\s*No\s+Component\s+Rate\s+Weight\s+Parameters", line):
            capture = True
            continue
        if capture:
            if not line.strip() or not re.search(r'\d', line):
                break
            parts = line.strip().split()
            if len(parts) >= 4:
                component = parts[1]
                weight = parts[3]
                profile_match = re.search(r'F(.+)', component)
                if profile_match:
                    profile = profile_match.group(1)
                    fmix_entries.append(f"{profile}:1.0:{weight}")

    return "FMIX{" + ",".join(fmix_entries) + "}" if fmix_entries else None

def main():
    if len(sys.argv) != 2:
        print("Usage: python extract_profiles_weights.py <input_directory>")
        sys.exit(1)

    input_dir = sys.argv[1]

    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    iqtree_files = [f for f in os.listdir(input_dir) if f.endswith(".iqtree")]
    output_files = []

    for iqfile in iqtree_files:
        full_path = os.path.join(input_dir, iqfile)
        fmix_result = extract_fmix_from_file(full_path)
        if fmix_result:
            prefix = os.path.splitext(iqfile)[0]
            out_path = os.path.join(input_dir, f"{prefix}.txt")
            with open(out_path, 'w', encoding='utf-8') as f_out:
                f_out.write(fmix_result + "\n")
            output_files.append(out_path)

    print(f"Done. {len(output_files)} FMIX profiles extracted.")

if __name__ == "__main__":
    main()

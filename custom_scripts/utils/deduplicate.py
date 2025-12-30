import argparse


def deduplicate_file(input_path, output_path):
    """Remove duplicate lines from a text file while preserving order."""
    seen = set()
    unique_lines = []
    
    with open(input_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped and stripped not in seen:
                seen.add(stripped)
                unique_lines.append(line)
    
    with open(output_path, "w") as f:
        f.writelines(unique_lines)
    
    total = len(open(input_path).readlines())
    print(f"Total lines: {total}")
    print(f"Unique lines: {len(unique_lines)}")
    print(f"Duplicates removed: {total - len(unique_lines)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove duplicate lines from a text file.")
    parser.add_argument("--input", "-i", required=True, help="Path to input file.")
    parser.add_argument("--output", "-o", required=True, help="Path to output file.")
    
    args = parser.parse_args()
    deduplicate_file(args.input, args.output)

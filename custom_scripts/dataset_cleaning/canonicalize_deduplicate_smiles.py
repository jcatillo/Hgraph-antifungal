import argparse
from rdkit import Chem


def canonicalize_smiles(input_path, output_path):
    """Convert SMILES to canonical form and drop duplicates."""
    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    total = len(smiles_list)
    seen = set()
    canonical = []
    invalid_count = 0
    duplicate_count = 0

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid_count += 1
            continue

        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        if canonical_smiles in seen:
            duplicate_count += 1
            continue

        seen.add(canonical_smiles)
        canonical.append(canonical_smiles)

    with open(output_path, "w") as f:
        for smiles in canonical:
            f.write(smiles + "\n")

    print("Completed canonicalization")
    print(f"  Total input: {total}")
    print(f"  Written (unique canonical): {len(canonical)}")
    print(f"  Dropped invalid: {invalid_count}")
    print(f"  Dropped duplicates after canonicalization: {duplicate_count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert SMILES to canonical SMILES and drop duplicates."
    )
    parser.add_argument("--input", "-i", required=True, help="Path to input SMILES file.")
    parser.add_argument(
        "--output", "-o", required=True, help="Path to write canonical unique SMILES."
    )

    args = parser.parse_args()
    canonicalize_smiles(args.input, args.output)

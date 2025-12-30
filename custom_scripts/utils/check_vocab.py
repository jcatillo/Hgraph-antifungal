import argparse
from rdkit import Chem


def check_vocab(filename):
    with open(filename, "r") as f:
        for i, line in enumerate(f, 1):
            parts = line.strip().split()
            if len(parts) < 2:
                print(f"Line {i}: Invalid format (less than 2 parts): {line.strip()}")
                continue

            smiles1, smiles2 = parts[0], parts[1]

            mol1 = Chem.MolFromSmiles(smiles1)
            if mol1 is None:
                print(f"Line {i}: Invalid SMILES 1: {smiles1}")

            mol2 = Chem.MolFromSmiles(smiles2)
            if mol2 is None:
                print(f"Line {i}: Invalid SMILES 2: {smiles2}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate SMILES pairs in a vocab file.")
    parser.add_argument("--input", "-i", required=True, help="Path to vocab file to validate.")

    args = parser.parse_args()
    check_vocab(args.input)

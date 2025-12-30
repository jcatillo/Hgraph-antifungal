import argparse
from rdkit import Chem


def analyze_dataset(input_path, output_path):
    with open(input_path, "r") as f:
        molecules = [line.strip() for line in f if line.strip()]

    total = len(molecules)
    print(f"Total molecules: {total}")

    disconnected_compounds = []
    clean_compounds = []
    invalid_smiles = []

    for smiles in molecules:
        if "." in smiles:
            disconnected_compounds.append(smiles)
        else:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                invalid_smiles.append(smiles)
            else:
                clean_compounds.append(smiles)

    print("\nBreakdown:")
    print(f"  Disconnected Compounds: {len(disconnected_compounds)}")
    print(f"  Connected valid: {len(clean_compounds)}")
    print(f"  Invalid SMILES: {len(invalid_smiles)}")

    pct_connected = (len(clean_compounds) / total * 100) if total else 0.0
    print(f"\nPercentage of connected compounds: {pct_connected:.1f}%")

    with open(output_path, "w") as f:
        for smiles in clean_compounds:
            f.write(smiles + "\n")

    print(f"\nSaved {len(clean_compounds)} connected compounds to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze dataset and save connected SMILES.")
    parser.add_argument("--input", "-i", required=True, help="Path to input SMILES file.")
    parser.add_argument("--output", "-o", required=True, help="Path to write connected SMILES.")

    args = parser.parse_args()
    analyze_dataset(args.input, args.output)
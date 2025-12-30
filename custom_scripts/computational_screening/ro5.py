import csv
import argparse
import os
from rdkit import Chem
from rdkit.Chem import Descriptors

# -------------------------------
# Load molecules
# -------------------------------
def load_smiles_from_csv(filename):
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        return [row for row in reader]

# -------------------------------
# Lipinski Rule of 5 check
# -------------------------------
def lipinski_check(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    # Count violations
    violations = 0
    violations += 1 if mw > 500 else 0
    violations += 1 if logp > 5 else 0
    violations += 1 if hbd > 5 else 0
    violations += 1 if hba > 10 else 0

    is_drug_like = violations <= 1  # True if 0 or 1 violation
    return mw, logp, hbd, hba, violations, is_drug_like


def screen_lipinski(input_csv, output_dir):
    """Screen molecules for Lipinski Rule of 5 compliance."""
    os.makedirs(output_dir, exist_ok=True)
    
    novel_mols = load_smiles_from_csv(input_csv)
    
    # -------------------------------
    # Process molecules
    # -------------------------------
    all_results = []
    drug_like_results = []

    for row in novel_mols:
        smi = row['canonical_smiles']
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            entry = [row['original_smiles'], smi, "INVALID", "INVALID", "INVALID", "INVALID", "INVALID", "INVALID"]
        else:
            mw, logp, hbd, hba, violations, is_drug_like = lipinski_check(mol)
            entry = [
                row['original_smiles'],
                smi,
                round(mw, 2),
                round(logp, 2),
                hbd,
                hba,
                violations,
                is_drug_like
            ]
        all_results.append(entry)
        if mol and is_drug_like:
            drug_like_results.append(entry)

    # -------------------------------
    # Save CSVs
    # -------------------------------
    headers = ["original_smiles", "canonical_smiles", "MW", "LogP", "HBD", "HBA", "violations", "is_drug_like"]
    
    ro5_all_csv = os.path.join(output_dir, "lipinski_all.csv")
    ro5_drug_like_csv = os.path.join(output_dir, "lipinski_drug_like.csv")
    ro5_drug_like_txt = os.path.join(output_dir, "drug_like_canonical_smiles.txt")

    with open(ro5_all_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(all_results)

    with open(ro5_drug_like_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(drug_like_results)

    # Save canonical SMILES of drug-like molecules to text file
    with open(ro5_drug_like_txt, "w") as f:
        for row in drug_like_results:
            canonical_smi = row[1]
            if canonical_smi != "INVALID":
                f.write(canonical_smi + "\n")

    print(f"Done! Saved all Lipinski results: {ro5_all_csv}")
    print(f"Done! Saved drug-like molecules CSV: {ro5_drug_like_csv}")
    print(f"Done! Saved drug-like canonical SMILES: {ro5_drug_like_txt}")
    print(f"Total screened: {len(all_results)}")
    print(f"Drug-like molecules (≤1 violation): {len(drug_like_results)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Screen molecules for Lipinski Rule of 5 compliance."
    )
    parser.add_argument("--input", "-i", required=True, help="Path to input CSV with canonical_smiles column.")
    parser.add_argument("--output", "-o", required=True, help="Output directory for Lipinski screening results.")

    args = parser.parse_args()
    screen_lipinski(args.input, args.output)
import csv
import argparse
import os
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# -------------------------------
# Load SMILES
# -------------------------------
def load_smiles(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip()]

# -------------------------------
# RDKit conversion helper
# -------------------------------
def to_mol(smi):
    mol = Chem.MolFromSmiles(smi)
    return mol if mol else None


def compute_tanimoto(generated_file, training_file, output_dir, threshold=0.7):
    """Compute Tanimoto similarity between generated and training molecules."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    gen_smiles = load_smiles(generated_file)
    train_smiles = load_smiles(training_file)

    train_fps = [AllChem.GetMorganFingerprintAsBitVect(to_mol(s), 2, nBits=2048)
                 for s in train_smiles if to_mol(s)]

    # -------------------------------
    # Process generated molecules
    # -------------------------------
    all_results = []
    novel_results = []

    for smi in gen_smiles:
        mol = to_mol(smi)
        if not mol:
            all_results.append([smi, "INVALID", "INVALID", "INVALID"])
            continue

        canon = Chem.MolToSmiles(mol, canonical=True)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

        # Compute Tanimoto similarities
        scores = [DataStructs.TanimotoSimilarity(fp, tfp) for tfp in train_fps]

        if scores:
            max_score = max(scores)
            max_index = scores.index(max_score)
            closest_smiles = train_smiles[max_index]
        else:
            max_score = 0.0
            closest_smiles = "NONE"

        row = [smi, canon, round(max_score, 4), closest_smiles]
        all_results.append(row)

        # Check novelty
        if max_score <= threshold:
            novel_results.append(row)

    # -------------------------------
    # Save CSVs
    # -------------------------------
    all_csv = os.path.join(output_dir, "tanimoto_all_gen.csv")
    novel_csv = os.path.join(output_dir, f"tanimoto_threshold_{threshold}.csv")
    novel_txt = os.path.join(output_dir, f"novel_canonical_smiles_{threshold}.txt")

    with open(all_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["original_smiles", "canonical_smiles", "tanimoto_score", "closest_training_molecule"])
        writer.writerows(all_results)

    with open(novel_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["original_smiles", "canonical_smiles", "tanimoto_score", "closest_training_molecule"])
        writer.writerows(novel_results)

    # Save canonical SMILES of novel molecules to text file
    with open(novel_txt, "w") as f:
        for row in novel_results:
            canonical_smi = row[1]
            if canonical_smi != "INVALID":
                f.write(canonical_smi + "\n")

    print(f"Done! Saved all molecules: {all_csv}")
    print(f"Done! Saved novel molecules CSV: {novel_csv}")
    print(f"Done! Saved canonical SMILES of novel molecules: {novel_txt}")
    print(f"Total generated: {len(gen_smiles)}")
    print(f"Novel molecules (Tanimoto ≤ {threshold}): {len(novel_results)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute Tanimoto similarity between generated and training molecules."
    )
    parser.add_argument("--input", "-i", required=True, help="Path to generated molecules file (SMILES, one per line).")
    parser.add_argument("--training", "-t", required=True, help="Path to training/dataset molecules file (SMILES, one per line).")
    parser.add_argument("--output", "-o", required=True, help="Output directory for result CSVs.")
    parser.add_argument("--threshold", type=float, default=0.7, help="Max Tanimoto similarity for novelty (default: 0.7).")

    args = parser.parse_args()
    compute_tanimoto(args.input, args.training, args.output, args.threshold)
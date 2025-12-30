import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import os

def calculate_properties(smiles):
    """Calculate Ro5 properties for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    
    # Calculate Violations
    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    
    return {
        'MW': mw,
        'LogP': logp,
        'HBD': hbd,
        'HBA': hba,
        'Violations': violations
    }

import argparse

def main():
    parser = argparse.ArgumentParser(description='Calculate physicochemical properties for molecules.')
    parser.add_argument('--input', '-i', default="antifungal_screening_multi_target_results.csv", help='Input CSV file path')
    parser.add_argument('--output', '-o', default="properties_output.csv", help='Output CSV file path')
    args = parser.parse_args()
    
    csv_file = args.input
    output_file = args.output
    
    if not os.path.exists(csv_file):
        print(f"Error: File {csv_file} not found.")
        return
        
    print(f"Loading data from {csv_file}...")
    df = pd.read_csv(csv_file)
    
    # Get all molecules sorted by Rank
    all_molecules = df.sort_values('Rank')
    
    results = []
    
    for _, row in all_molecules.iterrows():
        rank = row['Rank']
        smiles = row['SMILES']
        props = calculate_properties(smiles)
        
        if props:
            results.append({
                'Rank': rank,
                'SMILES': smiles,
                'MW': round(props['MW'], 2),
                'LogP': round(props['LogP'], 2),
                'HBD': props['HBD'],
                'HBA': props['HBA'],
                'Violations': props['Violations']
            })
    
    # Save to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")
    print(f"Total molecules processed: {len(results_df)}")

if __name__ == "__main__":
    main()

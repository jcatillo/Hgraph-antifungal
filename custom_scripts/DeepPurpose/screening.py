"""
Multi-Target DeepPurpose Screening
Screen molecules against multiple protein targets simultaneously
"""

import pandas as pd
import numpy as np
import os
from DeepPurpose import DTI as models
from DeepPurpose import utils
import argparse
from tqdm import tqdm


def load_targets(targets_file):
    """
    Load protein targets from file.
    Supports:
    1. CSV file with 'Entry Name' (or 'Entry') and 'Sequence' columns
    2. Text file with "TargetName: SEQUENCE" format
    3. Text file with just sequences
    """
    targets = {}
    
    if targets_file.endswith('.csv'):
        try:
            df = pd.read_csv(targets_file)
            # Check for required columns
            if 'Sequence' not in df.columns:
                raise ValueError("CSV must contain a 'Sequence' column")
            
            # Determine name column
            name_col = 'Entry Name' if 'Entry Name' in df.columns else 'Entry'
            if name_col not in df.columns:
                name_col = df.columns[0]  # Fallback to first column
            
            for _, row in df.iterrows():
                name = str(row[name_col]).strip()
                seq = str(row['Sequence']).strip()
                if seq and seq.lower() != 'nan':
                    targets[name] = seq
            
            print(f"Loaded {len(targets)} targets from CSV")
            return targets
        except Exception as e:
            print(f"Error reading CSV: {e}")
            print("Falling back to text parsing...")
            
    # Fallback to text parsing
    with open(targets_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    for i, line in enumerate(lines, 1):
        if ':' in line:
            name, seq = line.split(':', 1)
            targets[name.strip()] = seq.strip()
        else:
            targets[f'Target_{i}'] = line
    
    return targets


def multi_target_screening(smiles_file, targets_file, output_prefix, 
                           model_name='MPNN_CNN_BindingDB'):
    """
    Screen molecules against multiple protein targets.
    """
    print("=" * 80)
    print("Multi-Target Drug-Target Interaction Screening")
    print("=" * 80)
    print()
    
    # Load molecules
    print(f"Loading molecules from: {smiles_file}")
    with open(smiles_file, 'r') as f:
        smiles_list = [line.strip() for line in f if line.strip()]
    print(f"Total molecules: {len(smiles_list)}")
    print()
    
    # Load targets
    print(f"Loading protein targets from: {targets_file}")
    targets = load_targets(targets_file)
    print(f"Total targets: {len(targets)}")
    for name, seq in targets.items():
        print(f"  - {name}: {len(seq)} amino acids")
    print()
    
    # Load model once
    print(f"Loading pre-trained model: {model_name}...")
    # Set model directory to avoid creating save_folder in root
    model_dir = 'custom_scripts/DeepPurpose/pretrained_models'
    os.makedirs(model_dir, exist_ok=True)
    
    # Change to model directory temporarily
    original_dir = os.getcwd()
    os.chdir(model_dir)
    try:
        model = models.model_pretrained(model=model_name)
    finally:
        os.chdir(original_dir)
    print()
    
    # Results storage
    all_results = []
    
    # Screen against each target
    for target_name, target_seq in targets.items():
        print(f"Screening against: {target_name}")
        print("-" * 80)
        
        # Prepare data
        drug_list = smiles_list
        target_list = [target_seq] * len(smiles_list)
        
        # Encode
        X_pred = utils.data_process(
            X_drug=drug_list,
            X_target=target_list,
            y=[0] * len(drug_list),
            drug_encoding='MPNN',
            target_encoding='CNN',
            split_method='no_split'
        )
        
        # Predict
        predictions = model.predict(X_pred)
        
        # Store results
        for smiles, affinity in zip(smiles_list, predictions):
            all_results.append({
                'SMILES': smiles,
                'Target': target_name,
                'Predicted_Affinity': affinity
            })
        
        print(f"✓ Completed {target_name}")
        print(f"  Mean affinity: {np.mean(predictions):.3f}")
        print(f"  Max affinity: {np.max(predictions):.3f}")
        print()
    
    # Create results dataframe
    df = pd.DataFrame(all_results)
    
    # Pivot table: rows = molecules, columns = targets
    pivot_df = df.pivot(index='SMILES', columns='Target', values='Predicted_Affinity')
    
    # Add summary statistics
    pivot_df['Mean_Affinity'] = pivot_df.mean(axis=1)
    pivot_df['Max_Affinity'] = pivot_df.max(axis=1)
    pivot_df['Min_Affinity'] = pivot_df.min(axis=1)
    pivot_df['Std_Affinity'] = pivot_df.std(axis=1)
    
    # Sort by mean affinity
    pivot_df = pivot_df.sort_values('Mean_Affinity', ascending=False)
    pivot_df['Rank'] = range(1, len(pivot_df) + 1)
    
    # Reorder columns
    cols = ['Rank', 'Mean_Affinity', 'Max_Affinity', 'Min_Affinity', 'Std_Affinity'] + \
           [col for col in pivot_df.columns if col.startswith('Target_') or col in targets.keys()]
    pivot_df = pivot_df[[col for col in cols if col in pivot_df.columns]]
    
    # Summary
    print("=" * 80)
    print("MULTI-TARGET SCREENING RESULTS")
    print("=" * 80)
    print()
    print(f"Molecules screened: {len(smiles_list)}")
    print(f"Targets screened: {len(targets)}")
    print(f"Total predictions: {len(all_results)}")
    print()
    
    # Top 10 broad-spectrum candidates (high mean affinity)
    print("Top 10 Broad-Spectrum Candidates (by Mean Affinity):")
    print("-" * 80)
    for idx, row in pivot_df.head(10).iterrows():
        print(f"Rank {int(row['Rank']):2d}: Mean={row['Mean_Affinity']:.3f} | {idx[:60]}...")
    print()
    
    # Best binder for each target
    print("Best Binder for Each Target:")
    print("-" * 80)
    for target_name in targets.keys():
        if target_name in pivot_df.columns:
            best_idx = pivot_df[target_name].idxmax()
            best_affinity = pivot_df.loc[best_idx, target_name]
            print(f"{target_name}: {best_affinity:.3f} | {best_idx[:50]}...")
    print()
    
    # Save results
    output_file = f"{output_prefix}_multi_target_results.csv"
    pivot_df.to_csv(output_file)
    print(f"✓ Results saved to: {output_file}")
    
    # Save top 10 broad-spectrum
    top10_file = f"{output_prefix}_top10_broad_spectrum.txt"
    with open(top10_file, 'w') as f:
        for smiles in pivot_df.head(10).index:
            f.write(smiles + '\n')
    print(f"✓ Top 10 broad-spectrum SMILES saved to: {top10_file}")
    
    # Save detailed results (long format)
    detailed_file = f"{output_prefix}_detailed_results.csv"
    df.to_csv(detailed_file, index=False)
    print(f"✓ Detailed results saved to: {detailed_file}")
    
    return pivot_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Screen molecules against multiple protein targets'
    )
    parser.add_argument('--smiles', '-s', required=True, help='File with SMILES strings')
    parser.add_argument('--targets', '-t', required=True, help='File with protein sequences (CSV or txt)')
    parser.add_argument('--output', '-o', default='results/DeepPurpose/screening',
                       help='Output file prefix (default: results/dti_screening/screening)')
    parser.add_argument('--model', default='MPNN_CNN_BindingDB',
                       help='Pre-trained model (default: MPNN_CNN_BindingDB)')
    
    args = parser.parse_args()
    
    multi_target_screening(args.smiles, args.targets, 
                          args.output, args.model)
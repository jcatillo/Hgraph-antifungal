import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw
import os
import numpy as np

# Set style for plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("viridis")

def load_data(filepath):
    """Load the screening results CSV."""
    print(f"Loading data from {filepath}...")
    df = pd.read_csv(filepath)
    return df

def generate_molecule_images(df, output_dir="screening_molecules"):
    """Generate PNG images for each molecule."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"Generating molecule images in {output_dir}...")
    
    for index, row in df.iterrows():
        smiles = row['SMILES']
        rank = row['Rank']
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Create a valid filename
            filename = f"rank_{rank}.png"
            filepath = os.path.join(output_dir, filename)
            
            # Draw with caption (legend)
            # Truncate SMILES if too long for the image width, or let RDKit handle it (it might wrap or shrink)
            # For very long SMILES, RDKit's default font size might be too small or text too crowded.
            # But we will try the standard legend first.
            try:
                img = Draw.MolToImage(mol, size=(400, 400), legend=smiles)
                img.save(filepath)
            except Exception as e:
                print(f"Error drawing {rank}: {e}")
                # Fallback without legend if it fails
                Draw.MolToFile(mol, filepath, size=(300, 300))
        else:
            print(f"Warning: Could not parse SMILES for Rank {rank}")

def plot_ranking(df, output_file="screening_ranking.png"):
    """Create a visual ranking based on Mean Affinity."""
    print(f"Creating ranking plot: {output_file}...")
    
    # Sort by Rank just to be sure
    df_sorted = df.sort_values('Rank')
    
    plt.figure(figsize=(12, 8))
    # Create bar plot
    # Plot all molecules as requested
    n_to_plot = len(df)
    plot_data = df_sorted.head(n_to_plot)
    
    bars = plt.bar(plot_data['Rank'].astype(str), plot_data['Mean_Affinity'], 
                   yerr=plot_data['Std_Affinity'], capsize=5, 
                   color=sns.color_palette("viridis", n_to_plot))
    
    # Add SD text annotations
    for bar, sd in zip(bars, plot_data['Std_Affinity']):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + sd + 0.1,
                 f'±{sd:.2f}',
                 ha='center', va='bottom', fontsize=6, rotation=90, color='black')
    
    plt.xlabel('Rank', fontsize=12)
    plt.ylabel('Mean Affinity (pKd)', fontsize=12)
    plt.title(f'Ranking of All {n_to_plot} Candidates by Mean Affinity', fontsize=14)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_consistency(df, output_file="top3_consistency.png"):
    """Visualize consistency of top 3 candidates across all targets."""
    print(f"Creating consistency plot: {output_file}...")
    
    # Get top 3 candidates
    top3 = df.sort_values('Rank').head(3)
    
    # Identify target columns (those starting with 'A0A')
    target_cols = [col for col in df.columns if col.startswith('A0A')]
    
    if not target_cols:
        print("Error: No target columns found (starting with 'A0A')")
        return

    # Prepare data for plotting (melt/unpivot)
    plot_data = []
    
    for _, row in top3.iterrows():
        rank = row['Rank']
        smiles = row['SMILES']
        # Shorten SMILES for legend if needed, or just use Rank
        label = f"Rank {rank}"
        
        for target in target_cols:
            affinity = row[target]
            plot_data.append({
                'Candidate': label,
                'Target': target,
                'Affinity': affinity
            })
            
    plot_df = pd.DataFrame(plot_data)
    
    plt.figure(figsize=(14, 8))
    
    # Violin plot to show distribution density
    sns.violinplot(x='Candidate', y='Affinity', data=plot_df, inner='quartile', palette="muted")
    
    # Overlay strip plot to show individual points
    sns.stripplot(x='Candidate', y='Affinity', data=plot_df, color='black', alpha=0.3, jitter=True)
    
    plt.ylabel('Binding Affinity (pKd)', fontsize=12)
    plt.title('Consistency of Top 3 Candidates Across All Targets', fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add mean lines or stats if needed
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

import argparse

def main():
    parser = argparse.ArgumentParser(description='Visualize multi-target screening results.')
    parser.add_argument('--input', '-i', default="antifungal_screening_multi_target_results.csv", help='Input CSV file path')
    parser.add_argument('--output-dir', '-o', default="screening_molecules", help='Directory to save molecule images')
    parser.add_argument('--ranking-plot', default="screening_ranking.png", help='Filename for the ranking plot')
    parser.add_argument('--consistency-plot', default="top3_consistency.png", help='Filename for the consistency plot')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: File {args.input} not found.")
        return
        
    df = load_data(args.input)
    
    # 1. Generate Images
    generate_molecule_images(df, output_dir=args.output_dir)
    
    # 2. Visual Ranking
    plot_ranking(df, output_file=args.ranking_plot)
    
    # 3. Consistency Analysis
    plot_consistency(df, output_file=args.consistency_plot)
    
    print("\nVisualization complete!")

if __name__ == "__main__":
    main()

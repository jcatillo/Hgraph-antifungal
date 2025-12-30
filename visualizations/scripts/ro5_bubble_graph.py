#!/usr/bin/env python3
"""
Lipinski's Rule of 5 - Bubble/Scatter Plot with Multiple Encodings
X-axis: MW | Y-axis: LogP | Size: HBD | Color: HBA
Edge Color: Number of Violations
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import numpy as np
import argparse


def calculate_ro5_properties(smiles):
    """Calculate Lipinski's Rule of 5 properties."""
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    
    violations = 0
    if mw > 500:
        violations += 1
    if logp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1
    
    return {
        'MW': mw,
        'LogP': logp,
        'HBD': hbd,
        'HBA': hba,
        'Violations': violations,
        'Pass_Ro5': violations == 0
    }


def create_bubble_plot(csv_file, output_file='ro5_bubble.png', smiles_column='canonical_smiles'):
    """
    Create bubble/scatter plot with multiple encodings.
    X: MW, Y: LogP, Size: HBD, Color: HBA, Edge Color: Violations
    """
    
    # Read and process data
    df = pd.read_csv(csv_file)
    
    if smiles_column not in df.columns:
        print(f"Error: Column '{smiles_column}' not found")
        return
    
    print(f"Processing {len(df)} molecules...")
    
    results = []
    for idx, row in df.iterrows():
        props = calculate_ro5_properties(row[smiles_column])
        if props:
            props['SMILES'] = row[smiles_column]
            props['Index'] = idx
            results.append(props)
    
    ro5_df = pd.DataFrame(results)
    
    # Statistics
    total = len(ro5_df)
    acceptable = (ro5_df['Violations'] <= 1).sum()
    poor = (ro5_df['Violations'] >= 2).sum()
    
    print(f"\nTotal: {total} | Acceptable (≤1 viol): {acceptable} | Poor (≥2 viol): {poor}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Define edge colors based on violations
    violation_colors = {
        0: 'green',
        1: 'green',
        2: 'black',
        3: 'orange',  # Using orange for 3 violations
        4: 'red'
    }
    
    # Plot molecules grouped by violations
    for viol in sorted(ro5_df['Violations'].unique()):
        viol_df = ro5_df[ro5_df['Violations'] == viol]
        sizes = viol_df['HBD'] * 80 + 120
        edge_color = violation_colors.get(viol, 'gray')
        edge_width = 3 if viol <= 1 else 2
        
        scatter = ax.scatter(viol_df['MW'], viol_df['LogP'],
                           s=sizes,
                           c=viol_df['HBA'],
                           cmap='viridis',
                           alpha=0.6,
                           edgecolors=edge_color,
                           linewidth=edge_width,
                           vmin=0,
                           vmax=20,
                           marker='o')
        
        # Add checkmarks for acceptable molecules (≤1 violation)
        if viol <= 1:
            for idx, row in viol_df.iterrows():
                ax.text(row['MW'], row['LogP'], '✓', 
                       fontsize=8, color='green', fontweight='bold',
                       ha='center', va='center', zorder=10)
    
    # Threshold lines
    ax.axvline(x=500, color='red', linestyle='--', linewidth=3, 
              alpha=0.8, label='MW threshold (500 Da)', zorder=5)
    ax.axhline(y=5, color='blue', linestyle='--', linewidth=3, 
              alpha=0.8, label='LogP threshold (5)', zorder=5)
    
    # Shade drug-like region (MW ≤ 500 AND LogP ≤ 5)
    drug_like_box = Rectangle((0, ax.get_ylim()[0]), 500, 5 - ax.get_ylim()[0],
                             linewidth=0, edgecolor=None, 
                             facecolor='green', alpha=0.08, zorder=0)
    ax.add_patch(drug_like_box)
    
    # Labels and title
    ax.set_xlabel('Molecular Weight (Da)', fontsize=14, fontweight='bold')
    ax.set_ylabel('LogP (Lipophilicity)', fontsize=14, fontweight='bold')
    ax.set_title('Lipinski\'s Rule of 5 - Comprehensive Bubble Plot\nBubble Size = H-Bond Donors | Color = H-Bond Acceptors', 
                fontsize=16, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    plt.tight_layout(rect=[0, 0.22, 0.85, 0.96])  # Make room on right for legend
    
    # Create custom legend with violation color coding - OUTSIDE plot area
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markeredgecolor='green', markeredgewidth=3, markersize=10, 
               label='Acceptable (≤1 viol) ✓'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markeredgecolor='black', markeredgewidth=2, markersize=10, 
               label='2 violations'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markeredgecolor='orange', markeredgewidth=2, markersize=10, 
               label='3 violations'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markeredgecolor='red', markeredgewidth=2, markersize=10, 
               label='4 violations'),
        Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='MW threshold (500 Da)'),
        Line2D([0], [0], color='blue', linestyle='--', linewidth=2, label='LogP threshold (5)')
    ]
    
    # Place legend outside to the right (top position)
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1.0),
             fontsize=11, framealpha=0.95, title='Legend', title_fontsize=12)
    
    # Colorbar for HBA - OUTSIDE plot area (right side)
    cbar_ax = fig.add_axes([0.88, 0.28, 0.025, 0.4])  # [left, bottom, width, height]
    cbar = plt.colorbar(scatter, cax=cbar_ax)
    cbar.set_label('H-Bond Acceptors (HBA)', fontsize=11, fontweight='bold', rotation=270, labelpad=20)
    cbar.ax.tick_params(labelsize=9)
    
    # Size legend for HBD - OUTSIDE plot area with ACTUAL bubble sizes
    legend_ax = fig.add_axes([0.35, 0.13, 0.3, 0.06])  # [left, bottom, width, height]
    legend_ax.set_xlim(0, 5)
    legend_ax.set_ylim(0, 1)
    legend_ax.axis('off')
    
    # Draw actual sized bubbles
    hbd_values = [0, 2, 5, 8, 10]
    x_positions = [0.5, 1.25, 2, 2.75, 3.5]
    
    for i, hbd in enumerate(hbd_values):
        # Calculate actual size (matching the plot)
        size = hbd * 80 + 120
        # Convert to radius for display (matplotlib uses area, so sqrt)
        radius = np.sqrt(size) / 100  # Scale down for legend
        
        circle = plt.Circle((x_positions[i], 0.5), radius, 
                           color='gray', alpha=0.6, edgecolor='black', linewidth=1)
        legend_ax.add_patch(circle)
        legend_ax.text(x_positions[i], 0.1, f'{hbd}', 
                      ha='center', va='top', fontsize=9, fontweight='bold')
    
    legend_ax.text(2, 0.95, 'Bubble Size (H-Bond Donors)', 
                  ha='center', va='top', fontsize=11, fontweight='bold')
    
    # Stats box - OUTSIDE plot area (bottom position) - COMBINED STATS
    stats_text = f"""Dataset Summary: Total = {total} | Acceptable (≤1 viol) = {acceptable} ({acceptable/total*100:.1f}%) | Poor (≥2 viol) = {poor} ({poor/total*100:.1f}%)
Ro5 Criteria: MW ≤ 500 Da | LogP ≤ 5 | HBD ≤ 5 | HBA ≤ 10"""
    
    fig.text(0.5, 0.04, stats_text, ha='center', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', 
                     alpha=0.9, edgecolor='black', linewidth=1.5))
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nBubble plot saved to: {output_file}")
    
    return ro5_df


def main():
    parser = argparse.ArgumentParser(description='Bubble plot with multiple encodings for Ro5')
    parser.add_argument('--input', '-i', required=True, help='Input CSV with SMILES')
    parser.add_argument('--output', '-o', default='ro5_bubble.png', help='Output file')
    parser.add_argument('--smiles-column', '-s', default='canonical_smiles', help='SMILES column')
    
    args = parser.parse_args()
    create_bubble_plot(args.input, args.output, args.smiles_column)


if __name__ == '__main__':
    main()
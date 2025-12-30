#!/usr/bin/env python3
"""
Visualization script for Tanimoto similarity scores.
Creates a histogram with a threshold line to identify novel molecules.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse


def visualize_tanimoto_scores(csv_file, threshold=0.07, output_file=None):
    """
    Create a histogram of Tanimoto scores with a threshold line.
    
    Args:
        csv_file: Path to the CSV file containing Tanimoto scores
        threshold: Threshold value for novelty (default: 0.07)
        output_file: Optional path to save the figure
    """
    # Read the data
    df = pd.read_csv(csv_file)
    
    # Extract Tanimoto scores
    tanimoto_scores = df['tanimoto_score'].values
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Create histogram
    n, bins, patches = ax.hist(tanimoto_scores, bins=50, color='skyblue', 
                                edgecolor='black', alpha=0.7, label='Tanimoto Scores')
    
    # Add threshold line
    ax.axvline(x=threshold, color='red', linestyle='--', linewidth=2, 
               label=f'Novelty Threshold ({threshold})')
    
    # Calculate statistics
    total_molecules = len(tanimoto_scores)
    novel_molecules = np.sum(tanimoto_scores <= threshold)
    novel_percentage = (novel_molecules / total_molecules) * 100
    
    # Add statistics text box
    stats_text = f'Total Molecules: {total_molecules}\n'
    stats_text += f'Novel (≤{threshold}): {novel_molecules} ({novel_percentage:.1f}%)\n'
    stats_text += f'Known (>{threshold}): {total_molecules - novel_molecules} ({100-novel_percentage:.1f}%)'
    
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Customize plot
    ax.set_xlabel('Tanimoto Similarity Score', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency (Number of Molecules)', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Tanimoto Similarity Scores\n(Comparison with Training Dataset)', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Set x-axis limits
    ax.set_xlim(-0.05, 1.05)
    
    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {output_file}")
    else:
        plt.show()
    
    # Print summary statistics
    print("\n=== Tanimoto Score Analysis ===")
    print(f"Total molecules analyzed: {total_molecules}")
    print(f"Novel molecules (≤{threshold}): {novel_molecules} ({novel_percentage:.2f}%)")
    print(f"Known molecules (>{threshold}): {total_molecules - novel_molecules} ({100-novel_percentage:.2f}%)")
    print(f"Mean Tanimoto score: {np.mean(tanimoto_scores):.4f}")
    print(f"Median Tanimoto score: {np.median(tanimoto_scores):.4f}")
    print(f"Min Tanimoto score: {np.min(tanimoto_scores):.4f}")
    print(f"Max Tanimoto score: {np.max(tanimoto_scores):.4f}")
    print(f"Std deviation: {np.std(tanimoto_scores):.4f}")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize Tanimoto similarity scores with histogram and threshold line'
    )
    parser.add_argument('--input', '-i', type=str, default='tanimoto_all.csv',
                        help='Input CSV file with Tanimoto scores (default: tanimoto_all.csv)')
    parser.add_argument('--threshold', '-t', type=float, default=0.7,
                        help='Novelty threshold value (default: 0.7)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output file path for saving the figure (optional)')
    
    args = parser.parse_args()
    
    visualize_tanimoto_scores(args.input, args.threshold, args.output)


if __name__ == '__main__':
    main()

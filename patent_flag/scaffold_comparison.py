"""
Scaffold Comparison Visualization
Compares the number of non-empty entries in original data vs scaffold outputs
"""

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse

def count_non_empty_lines(file_path):
    """Count non-empty lines in a file"""
    if not Path(file_path).exists():
        return 0
    
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        count = sum(1 for line in f if line.strip())
    return count

def main():
    parser = argparse.ArgumentParser(
        description='Compare original data with scaffold extraction results and generate visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scaffold_comparison.py
  python scaffold_comparison.py --data data --scaffold scaffold --output results.png
  python scaffold_comparison.py -d ./data -s ./scaffold -o output/comparison.png
        """
    )
    
    parser.add_argument(
        '-d', '--data',
        type=str,
        default='data',
        help='Path to data directory (default: data)'
    )
    parser.add_argument(
        '-s', '--scaffold',
        type=str,
        default='scaffold',
        help='Path to scaffold directory (default: scaffold)'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='scaffold_comparison.png',
        help='Output path for visualization (default: scaffold_comparison.png)'
    )
    
    args = parser.parse_args()
    
    data_dir = Path(args.data)
    scaffold_dir = Path(args.scaffold)
    output_file = args.output
    
    # File pairs to compare
    files_info = {
        'generated_set': [
            ('tanimoto_ro5_compliant_molecules', 'tanimoto_ro5_compliant_molecules'),
            ('all_generated_mol', 'all_generated_mol'),
            ('top10_molecules', 'top10_molecules'),
        ],
        'training_dataset': [
            ('training_dataset', 'training_dataset'),
        ]
    }
    
    results = []
    
    print("\nComparing Data Files with Scaffolds")
    print("="*70)
    
    for category, file_pairs in files_info.items():
        for data_name, scaffold_name in file_pairs:
            data_file = data_dir / category / f"{data_name}.txt"
            scaffold_file = scaffold_dir / category / f"{scaffold_name}_scaffold.txt"
            
            data_count = count_non_empty_lines(data_file)
            scaffold_count = count_non_empty_lines(scaffold_file)
            
            results.append({
                'File': data_name.replace('_', ' ').title(),
                'Category': category.replace('_', ' ').title(),
                'Original Data': data_count,
                'Scaffolds Found': scaffold_count,
                'Not Found': data_count - scaffold_count
            })
            
            print(f"{data_name:40} | Data: {data_count:4d} | Scaffold: {scaffold_count:4d} | Missing: {data_count - scaffold_count:4d}")
    
    df = pd.DataFrame(results)
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Side-by-side comparison
    x = range(len(df))
    width = 0.35
    
    ax1 = axes[0]
    ax1.bar([i - width/2 for i in x], df['Original Data'], width, label='Original Data', color='steelblue', alpha=0.8)
    ax1.bar([i + width/2 for i in x], df['Scaffolds Found'], width, label='Scaffolds Found', color='coral', alpha=0.8)
    
    ax1.set_xlabel('Dataset', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Number of Entries', fontsize=11, fontweight='bold')
    ax1.set_title('Scaffold Extraction Comparison', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(df['File'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Stacked bar chart showing completion rate
    ax2 = axes[1]
    ax2.bar(x, df['Scaffolds Found'], label='Scaffolds Found', color='coral', alpha=0.8)
    ax2.bar(x, df['Not Found'], bottom=df['Scaffolds Found'], label='Not Found (Empty)', color='lightcoral', alpha=0.8)
    
    ax2.set_xlabel('Dataset', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Number of Entries', fontsize=11, fontweight='bold')
    ax2.set_title('Scaffold Extraction Breakdown', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(df['File'], rotation=45, ha='right')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    
    # Add percentage labels
    for i, row in df.iterrows():
        if row['Original Data'] > 0:
            percentage = (row['Scaffolds Found'] / row['Original Data']) * 100
            ax2.text(i, row['Original Data'] + 0.5, f'{percentage:.1f}%', 
                    ha='center', va='bottom', fontweight='bold', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Visualization saved as '{output_file}'")
    
    # Print summary statistics
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    total_original = df['Original Data'].sum()
    total_scaffolds = df['Scaffolds Found'].sum()
    total_missing = df['Not Found'].sum()
    overall_rate = (total_scaffolds / total_original * 100) if total_original > 0 else 0
    
    print(f"Total Original Entries: {total_original}")
    print(f"Total Scaffolds Found: {total_scaffolds}")
    print(f"Total Not Found: {total_missing}")
    print(f"Overall Success Rate: {overall_rate:.2f}%")
    print("="*70)
    
    return df

if __name__ == '__main__':
    main()

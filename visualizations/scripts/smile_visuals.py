#!/usr/bin/env python3
"""
Generate individual image files for each molecule comparison.
Creates separate files for generated and training molecules.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os
import argparse
from pathlib import Path


def draw_single_molecule(smiles, label, output_file, img_size=(600, 500)):
    """
    Draw a single molecule and save to file.
    
    Args:
        smiles: SMILES string
        label: Label text to display
        output_file: Output file path
        img_size: Image size (width, height)
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        print(f"Warning: Could not parse SMILES: {smiles}")
        return False
    
    # Create drawer
    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    
    # Draw molecule with label
    drawer.DrawMolecule(mol, legend=label)
    drawer.FinishDrawing()
    
    # Save to file
    with open(output_file, 'wb') as f:
        f.write(drawer.GetDrawingText())
    
    return True


def sanitize_filename(smiles, max_length=50):
    """
    Create a safe filename from SMILES string.
    """
    # Replace problematic characters
    safe = smiles.replace('/', '_').replace('\\', '_').replace(':', '_')
    safe = safe.replace('*', 'x').replace('?', '_').replace('"', '_')
    safe = safe.replace('<', '_').replace('>', '_').replace('|', '_')
    
    # Truncate if too long
    if len(safe) > max_length:
        safe = safe[:max_length]
    
    return safe


def generate_all_comparisons(csv_file, output_dir="molecule_images"):
    """
    Generate individual image files for all molecules in the CSV.
    
    Args:
        csv_file: Path to CSV file with columns: canonical_smiles, closest_training_molecule, tanimoto_score
        output_dir: Directory to save images
    """
    # Read data
    df = pd.read_csv(csv_file)
    print(f"Processing {len(df)} molecules from {csv_file}")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Create subdirectories
    generated_dir = output_path / "generated"
    training_dir = output_path / "training"
    generated_dir.mkdir(exist_ok=True)
    training_dir.mkdir(exist_ok=True)
    
    # Process each molecule
    success_count = 0
    failed = []
    
    for idx, row in df.iterrows():
        try:
            # Get data
            gen_smiles = row['canonical_smiles']
            train_smiles = row['closest_training_molecule']
            tanimoto = row['tanimoto_score']
            
            # Create base filename using index and truncated SMILES
            base_name = f"{idx:03d}_{sanitize_filename(gen_smiles, 30)}"
            
            # Labels
            gen_label = f"Generated\nSMILES: {gen_smiles}"
            train_label = f"Training (Closest Match)\nSMILES: {train_smiles}"
            
            # Output files
            gen_file = generated_dir / f"{base_name}.png"
            train_file = training_dir / f"{base_name}.png"
            
            # Draw molecules
            gen_success = draw_single_molecule(gen_smiles, gen_label, gen_file)
            train_success = draw_single_molecule(train_smiles, train_label, train_file)
            
            if gen_success and train_success:
                success_count += 1
                if (idx + 1) % 10 == 0:
                    print(f"Processed {idx + 1}/{len(df)} molecules...")
            else:
                failed.append((idx, gen_smiles))
                
        except Exception as e:
            print(f"Error processing molecule {idx}: {e}")
            failed.append((idx, str(e)))
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"{'='*60}")
    print(f"Total molecules: {len(df)}")
    print(f"Successfully processed: {success_count}")
    print(f"Failed: {len(failed)}")
    print(f"\nOutput directories:")
    print(f"  Generated molecules: {generated_dir.absolute()}")
    print(f"  Training molecules:  {training_dir.absolute()}")
    
    if failed:
        print(f"\nFailed molecules:")
        for idx, info in failed[:5]:  # Show first 5
            print(f"  Index {idx}: {info}")
        if len(failed) > 5:
            print(f"  ... and {len(failed) - 5} more")
    
    # Create index file
    index_file = output_path / "index.txt"
    with open(index_file, 'w') as f:
        f.write("Molecular Comparison Index\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total molecules: {len(df)}\n")
        f.write(f"Successfully processed: {success_count}\n\n")
        f.write("Index | Generated SMILES | Tanimoto | Files\n")
        f.write("-"*80 + "\n")
        
        for idx, row in df.iterrows():
            gen_smiles = row['canonical_smiles']
            tanimoto = row['tanimoto_score']
            base_name = f"{idx:03d}_{sanitize_filename(gen_smiles, 30)}"
            f.write(f"{idx:03d} | {gen_smiles[:40]:40s} | {tanimoto:.4f} | {base_name}.png\n")
    
    print(f"\nIndex file created: {index_file.absolute()}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate individual image files for molecule comparisons (generated vs closest training).'
    )
    parser.add_argument('--input', '-i', type=str, default='tanimoto_novel.csv',
                       help='Input CSV file (default: tanimoto_novel.csv)')
    parser.add_argument('--output', '-o', type=str, default='molecule_images',
                       help='Output directory (default: molecule_images)')
    
    args = parser.parse_args()
    
    generate_all_comparisons(args.input, args.output)


if __name__ == '__main__':
    main()

import argparse
from rdkit import Chem

def get_atom_types(filename):
    atom_types = set()
    with open(filename, 'r') as f:
        for line in f:
            smiles = line.strip().split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                for atom in mol.GetAtoms():
                    atom_types.add((atom.GetSymbol(), atom.GetFormalCharge()))
    return atom_types

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List unique atom types in a SMILES file.")
    parser.add_argument("--input", "-i", required=True, help="Path to input SMILES file.")

    args = parser.parse_args()
    atoms = get_atom_types(args.input)
    print("Found atoms:")
    for symbol, charge in sorted(atoms):
        print(f"('{symbol}', {charge})")

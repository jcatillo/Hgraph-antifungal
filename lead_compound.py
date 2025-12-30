from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# Your generated molecule SMILES
smiles = "C=CCCCCCCCCC(=O)O" 
mol = Chem.MolFromSmiles(smiles)

# Extract the scaffold
scaffold = MurckoScaffold.GetScaffoldForMol(mol)
scaffold_smiles = Chem.MolToSmiles(scaffold)

print(f"Original: {smiles}")
print(f"Scaffold: {scaffold_smiles}")
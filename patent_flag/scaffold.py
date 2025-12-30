import argparse
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def extract_scaffolds(input_txt, output_txt):
    try:
        with open(input_txt, 'r') as f:
            smiles_list = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: {input_txt} not found.")
        return

    lead_compounds = []

    for i, smi in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(smi)
            # If RDKit cannot parse the SMILES, it returns None
            if mol is None:
                lead_compounds.append("")
                continue

            # Get the Bemis-Murcko scaffold
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smi = Chem.MolToSmiles(scaffold) if scaffold else ""
            lead_compounds.append(scaffold_smi)

        except Exception:
            print(f"Row {i} failed")
            lead_compounds.append("")

    with open(output_txt, 'w') as f:
        for scaffold_smi in lead_compounds:
            f.write(scaffold_smi + "\n")

def main():
    parser = argparse.ArgumentParser(description="Generate Murcko scaffolds from exact SMILES")
    parser.add_argument("-i", "--input", required=True, help="Input TXT file")
    parser.add_argument("-o", "--output", required=True, help="Output TXT file")

    args = parser.parse_args()
    extract_scaffolds(args.input, args.output)

if __name__ == "__main__":
    main()
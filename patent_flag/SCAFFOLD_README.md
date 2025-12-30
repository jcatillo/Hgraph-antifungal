# Scaffold

Extracts Bemis-Murcko scaffolds from SMILES strings.

## Install

```bash
pip install rdkit
```

## Usage

```bash
python scaffold.py -i input.txt -o output.txt
```

### Arguments

- `-i, --input` (required): Input file with SMILES (one per line)
- `-o, --output` (required): Output file for scaffold SMILES

## Quick Examples

```bash
python scaffold.py -i molecules.txt -o scaffolds.txt
```

## Input Format

One SMILES per line:

```
CCO
c1ccccc1O
CC(C)Cc1ccc(cc1)C(C)C(O)=O
```

## Output Format

One scaffold SMILES per line (matching input order):

```
CO
c1ccccc1
c1ccc2ccccc2c1
```

## How It Works

1. Parses each input SMILES using RDKit
2. Extracts Bemis-Murcko scaffold (generic structure)
3. Converts scaffold to SMILES
4. Writes scaffolds to output file (preserves row count)
5. Skips invalid SMILES (outputs empty line)

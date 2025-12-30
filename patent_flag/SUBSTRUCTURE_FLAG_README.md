# Substructure Flag

Searches for substructure matches in PubChem and flags if any are patented.

## Install

```bash
pip install pandas requests
```

## Usage

```bash
python substructure_flag.py -i input.txt -o output.csv
```

### Arguments

- `-i, --input` (required): Input file with SMILES (one per line)
- `-o, --output` (default: `substructure_results.csv`): Output CSV file

## Quick Examples

```bash
python substructure_flag.py -i molecules.txt
python substructure_flag.py -i molecules.txt -o results.csv
```

## Input Format

One SMILES per line:

```
CCO
c1ccccc1O
CC(C)Cc1ccc(cc1)C
```

## Output

CSV with columns:

- `input_smile`: Input SMILES
- `is_patented`: Patent found (True/False)

## How It Works

1. Searches PubChem for compounds matching input as substructure
2. Uses async polling for results (handles long-running searches)
3. Checks top 50 results for patent information
4. Flags first patented match found
5. Respects API rate limits (0.2s delays, 2s polling intervals)

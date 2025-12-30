# Similar Flag

Flags molecular compounds (SMILES) by finding similar patented compounds in PubChem.

## Install

```bash
pip install pandas requests
```

## Usage

```bash
python similar_flag.py -i input.txt -o output.csv -t 95
```

### Arguments

- `-i, --input` (required): Input file with SMILES strings (one per line)
- `-o, --output` (default: `patent_results.csv`): Output CSV file
- `-t, --threshold` (default: `95`): Similarity threshold (0-100)

## Examples

```bash
python similar_flag.py -i molecules.txt
python similar_flag.py -i molecules.txt -o results.csv -t 85
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

1. Searches PubChem for 2D-similar compounds
2. Checks each result for patent information
3. Flags first patented match found
4. Respects API rate limits (0.2s delays)

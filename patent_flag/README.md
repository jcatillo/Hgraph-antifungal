# Patent Flag - Molecular Patent Analysis

Python toolkit for screening compounds against PubChem patents using scaffold extraction, similarity search, and substructure matching.

## Table of Contents

- [Quick Start](#quick-start)
- [Installation](#installation)
- [Scripts](#scripts)
- [Usage Examples](#usage-examples)
- [Output Format](#output-format)
- [File Structure](#file-structure)

---

## Quick Start

```bash
# 1. Extract scaffolds
python scaffold.py -i exact/molecules.txt -o scaffold/molecules.txt

# 2. Search similar compounds for patents
python similar_flag.py -i exact/molecules.txt -o similarity/results.csv -t 95

# 3. Search substructure matches for patents
python substructure_flag.py -i exact/molecules.txt -o substructure/results.csv
```

---

## Installation

```bash
pip install pandas requests rdkit matplotlib
```

---

## Scripts

| Script                 | Purpose                            | Command                                                  |
| ---------------------- | ---------------------------------- | -------------------------------------------------------- |
| `scaffold.py`          | Extract Bemis-Murcko scaffolds     | `python scaffold.py -i INPUT -o OUTPUT`                  |
| `similar_flag.py`      | Find 2D-similar patented compounds | `python similar_flag.py -i INPUT -o OUTPUT -t THRESHOLD` |
| `substructure_flag.py` | Find substructure-matching patents | `python substructure_flag.py -i INPUT -o OUTPUT`         |

---

## Usage Examples

### Scaffold Extraction (for Data input)

```bash
python scaffold.py -i data/generated_set/tanimoto_ro5_compliant_molecules.txt -o scaffold/generated_set/tanimoto_ro5_compliant_molecules_scaffold.txt
python scaffold.py -i data/generated_set/all_generated_mol.txt -o scaffold/generated_set/all_generated_mol_scaffold.txt
python scaffold.py -i data/generated_set/top10_molecules.txt -o scaffold/generated_set/top10_molecules_scaffold.txt
python scaffold.py -i data/training_dataset/training_dataset.txt -o scaffold/training_dataset/training_dataset_scaffold.txt
```

### Similarity Search - EXACT

```bash
python similar_flag.py -i data/generated_set/tanimoto_ro5_compliant_molecules.txt -o similarity/generated_set/tanimoto_ro5_compliant_molecules_exact_flag.csv -t 95
python similar_flag.py -i data/generated_set/all_generated_mol.txt -o similarity/generated_set/all_generated_mol_exact_flag.csv -t 95
python similar_flag.py -i data/generated_set/top10_molecules.txt -o similarity/generated_set/top_10_molecules_exact_flag.csv -t 95
python similar_flag.py -i data/training_dataset/training_dataset.txt -o similarity/training_dataset/training_dataset_exact_flag.csv -t 95
```

### Similarity Search - SCAFFOLD

```bash
python similar_flag.py -i scaffold/generated_set/tanimoto_ro5_compliant_molecules_scaffold.txt -o similarity/generated_set/tanimoto_ro5_compliant_molecules_scaffold_flag.csv -t 95
python similar_flag.py -i scaffold/generated_set/all_generated_mol_scaffold.txt -o similarity/generated_set/all_generated_mol_scaffold_flag.csv -t 95
python similar_flag.py -i scaffold/generated_set/top10_molecules_scaffold.txt -o similarity/generated_set/top10_molecules_scaffold_flag.csv -t 95
python similar_flag.py -i scaffold/training_dataset/training_dataset_scaffold.txt -o similarity/training_dataset/training_dataset_scaffold_flag.csv -t 95
```

### Substructure Search - EXACT

```bash
python substructure_flag.py -i data/generated_set/tanimoto_ro5_compliant_molecules.txt -o substructure/generated_set/tanimoto_ro5_compliant_molecules_exact_exact_flag.csv
python substructure_flag.py -i data/generated_set/all_generated_mol.txt -o substructure/generated_set/all_generated_mol_exact_flag.csv
python substructure_flag.py -i data/generated_set/top10_molecules.txt -o substructure/generated_set/top_10_molecules_exact_flag.csv
python substructure_flag.py -i data/training_dataset/training_dataset.txt -o substructure/training_dataset/training_dataset_exact_flag.csv
```

### Substructure Search - SCAFFOLD

```bash
python substructure_flag.py -i scaffold/generated_set/tanimoto_ro5_compliant_molecules_scaffold.txt -o substructure/generated_set/tanimoto_ro5_compliant_molecules_scaffold_flag.csv
python substructure_flag.py -i scaffold/generated_set/all_generated_mol_scaffold.txt -o substructure/generated_set/all_generated_mol_scaffold_flag.csv
python substructure_flag.py -i scaffold/generated_set/top10_molecules_scaffold.txt -o substructure/generated_set/top_10_molecules_scaffold_flag.csv
python substructure_flag.py -i scaffold/training_dataset/training_dataset_scaffold.txt -o substructure/training_dataset/training_dataset_scaffold_flag.csv
```

---

## Output Format

All CSV outputs contain:

- `input_smile`: Original SMILES string
- `is_patented`: True/False (patent found)

Example:

```csv
input_smile,is_patented
CCO,True
c1ccccc1O,False
```

---

## File Structure

```
patent_flag/
├── scaffold.py                          # Main scripts
├── similar_flag.py
├── substructure_flag.py
│
├── data/                                # Input SMILES files
│   ├── generated_set/
│   │   ├── tanimoto_ro5_compliant_molecules.txt
│   │   ├── all_generated_mol.txt
│   │   └── top10_molecules.txt
│   └── training_dataset/
│       └── training_dataset.txt
│
├── scaffold/                            # Generated scaffolds
│   ├── generated_set/
│   │   ├── tanimoto_ro5_compliant_molecules_scaffold.txt
│   │   ├── all_generated_mol_scaffold.txt
│   │   └── top10_molecules_scaffold.txt
│   └── training_dataset/
│       └── training_dataset_scaffold.txt
│
├── similarity/                          # Similarity search results
│   ├── generated_set/
│   │   ├── tanimoto_ro5_compliant_molecules_exact_flag.csv
│   │   ├── tanimoto_ro5_compliant_molecules_scaffold_flag.csv
│   │   ├── all_generated_mol_exact_flag.csv
│   │   ├── all_generated_mol_scaffold_flag.csv
│   │   ├── top_10_molecules_exact_flag.csv
│   │   └── top10_molecules_scaffold_flag.csv
│   └── training_dataset/
│       ├── training_dataset_exact_flag.csv
│       └── training_dataset_scaffold_flag.csv
│
└── substructure/                        # Substructure search results
    ├── generated_set/
    │   ├── tanimoto_ro5_compliant_molecules_exact_exact_flag.csv
    │   ├── tanimoto_ro5_compliant_molecules_scaffold_flag.csv
    │   ├── all_generated_mol_exact_flag.csv
    │   ├── all_generated_mol_scaffold_flag.csv
    │   ├── top_10_molecules_exact_flag.csv
    │   └── top_10_molecules_scaffold_flag.csv
    └── training_dataset/
        ├── training_dataset_exact_flag.csv
        └── training_dataset_scaffold_flag.csv
```

---

## Visualization Scripts

Generate comparison charts for analysis results:

### Scaffold Comparison

```bash
python visualizations/scaffold_comparison.py
python visualizations/scaffold_comparison.py -d data -s scaffold -o output.png
```

Compares non-empty entries in original data vs scaffold outputs.

### Similarity Search Comparison

```bash
python visualizations/similarity_comparison.py
python visualizations/similarity_comparison.py -s similarity -o output.png
```

Visualizes patent hit rates for exact vs scaffold similarity searches (Tanimoto threshold: 95%).

### Substructure Search Comparison

```bash
python visualizations/substructure_comparison.py
python visualizations/substructure_comparison.py -s substructure -o output.png
```

Visualizes patent hit rates for exact vs scaffold substructure searches.

---

**scaffold.py**

- `-i, --input` (required): Input TXT file
- `-o, --output` (required): Output TXT file

**similar_flag.py**

- `-i, --input` (required): Input TXT file
- `-o, --output` (default: `patent_results.csv`): Output CSV
- `-t, --threshold` (default: `95`): Similarity threshold 0-100

**substructure_flag.py**

- `-i, --input` (required): Input TXT file
- `-o, --output` (default: `substructure_results.csv`): Output CSV

---

## Input Format

One SMILES per line:

```
CCO
c1ccccc1O
CC(C)Cc1ccc(cc1)C(C)C(O)=O
```

---

## Notes

- Respects PubChem API rate limits (0.2s delays)
- Substructure search uses async polling (2s intervals)
- Empty lines or comments skipped automatically

---

## Documentation

- [SCAFFOLD_README.md](SCAFFOLD_README.md) - Scaffold tool details
- [SIMILAR_FLAG_README.md](SIMILAR_FLAG_README.md) - Similarity search details
- [SUBSTRUCTURE_FLAG_README.md](SUBSTRUCTURE_FLAG_README.md) - Substructure search details
- [GRAPH_INTERPRETATION.md](GRAPH_INTERPRETATION.md) - Analysis guide

---

**Last Updated**: December 2025

# Antifungal Pipeline Instructions

This document provides step-by-step instructions for running the antifungal workflow on hgraph2graph, with explicit references to the scripts and outputs used in each stage.

---

## Table of Contents
1. [Prereqs](#0-prereqs)
2. [Environment](#1-create-environment-python-310)
3. [Data Cleaning](#2-data-cleaning)
4. [Vocabulary](#3-vocabulary)
5. [Preprocess](#4-preprocess)
6. [Train Generator](#5-train-generator)
7. [Generate Molecules](#6-generate-molecules)
8. [Deduplicate & Canonicalize](#7-deduplicate--canonicalize-generated-set)
9. [Computational Screening](#8-computational-screening)
10. [CLI Parsers (Per Script)](#cli-parsers-per-script)
11. [Outputs Reference](#outputs-reference)
12. [Notes](#notes)

## 0) Prereqs
- Windows with PowerShell (or Linux/macOS shell); git installed.
- GPU with CUDA 12.8 for training/inference (optional but recommended). CPU works but slower.

## 1) Create Environment (Python 3.10)
```bash
conda create -n hgraph_env python=3.10 -y
conda activate hgraph_env
```

Install core deps (RDKit must come from conda-forge):
```bash
conda install -c conda-forge rdkit -y
pip install torch networkx tqdm numpy
pip install git+https://github.com/bp-kelley/descriptastorus
pip install DeepPurpose
```

If you need GPU PyTorch, install the CUDA 12.8 build from PyTorch instructions before the other pip deps.

## 2) Data Cleaning
Goal: produce single-component, canonical SMILES at data/antifungal/clean.txt.
- Canonicalize raw SMILES.
- Remove disconnected/fragmented molecules.
Use your preferred cleaner; example placeholder (replace with your script):
```bash
python custom_scripts/dataset_cleaning/canonicalize_smiles.py --input data/antifungal/all.txt --output data/antifungal/clean.txt
```

## 3) Vocabulary
Build vocab with 8 CPUs from the cleaned set:
```bash
python get_vocab.py --ncpu 8 < data/antifungal/clean.txt > data/antifungal/vocab.txt
```

## 4) Preprocess
Create training tensors from clean.txt:
```bash
python preprocess.py --train data/antifungal/clean.txt --vocab data/antifungal/vocab.txt --ncpu 8 --mode single
```
Outputs: tensor* files (place or keep in train_processed/ as needed).

## 5) Train Generator
```bash
python train_generator.py \
    --train test/ \
    --vocab data/antifungal/vocab.txt \
    --save_dir ckpt/antifungal \
    --batch_size 32 \
    --epoch 150 \
    --warmup 500 \
    --step_beta 0.005 \
    --kl_anneal_iter 10 \
    --anneal_iter 1000 \
    --print_iter 10 \
    --save_iter 100
```
Adjust --train to your preprocessed tensor directory if different.

## 6) Generate Molecules
```bash
python generate.py --vocab data/antifungal/vocab.txt --model ckpt/antifungal/model.ckpt.3300 --nsamples 1000
```
Output: gen.txt (SMILES).

## 7) Deduplicate & Canonicalize Generated Set
```bash
python custom_scripts/dataset_cleaning/canonicalize_smiles.py --input gen.txt --output gen_unique.txt
```

## 8) Computational Screening

### Novelty (Tanimoto <= 0.7)
```bash
python custom_scripts/computational_screening/tanimoto.py -i results/hgraph_gen/gen.txt -t data/antifungal/clean.txt -o results/computational_screening/tanimoto --threshold 0.7
```
Outputs:
- tanimoto_all_gen.csv
- tanimoto_threshold_0.7.csv (novel set)
- novel_canonical_smiles_0.7.txt

### Lipinski (RO5)
```bash
python custom_scripts/computational_screening/ro5.py -i results/computational_screening/tanimoto/tanimoto_threshold_0.7.csv -o results/computational_screening/lipinski
```
Outputs:
- lipinski_all.csv
- lipinski_drug_like.csv
- drug_like_canonical_smiles.txt

### DTI (DeepPurpose)
```bash
python custom_scripts/DeepPurpose/screening.py --smiles results/computational_screening/lipinski/drug_like_canonical_smiles.txt --targets data/antifungal/targets.csv --output results/dti_screening/screening
```
Outputs under results/dti_screening/: screening_multi_target_results.csv, screening_detailed_results.csv, screening_top10_broad_spectrum.txt

## Notes
- Paths are relative to repo root; run commands from hgraph2graph/.
- Ensure CUDA toolkit matches your GPU driver if using GPU builds.
- If DeepPurpose downloads models, they will cache under its default directory; adjust if needed.

## Outputs Reference
- Generation: `gen.txt` (root) from `generate.py`.
- Dedup/canonical: `gen_unique.txt` from `custom_scripts/dataset_cleaning/canonicalize_smiles.py`.
- Novelty: `results/computational_screening/tanimoto/`
    - `tanimoto_all_gen.csv`
    - `tanimoto_threshold_0.7.csv`
    - `novel_canonical_smiles_0.7.txt`
- Lipinski: `results/computational_screening/lipinski/`
    - `lipinski_all.csv`
    - `lipinski_drug_like.csv`
    - `drug_like_canonical_smiles.txt`
- DTI: `results/dti_screening/`
    - `screening_multi_target_results.csv`
    - `screening_detailed_results.csv`
    - `screening_top10_broad_spectrum.txt`

## CLI Parsers (Per Script)
Each script defines its own argparse interface. Required flags must be provided or the script exits.

- `get_vocab.py`
    - `--ncpu <int>`: number of CPU workers. Input via stdin; vocab.txt via stdout redirection.

- `preprocess.py`
    - `--train <file>`: training SMILES file.
    - `--vocab <file>`: vocab file.
    - `--ncpu <int>`: CPU workers.
    - `--mode <single|pair>`: data mode.

- `train_generator.py`
    - `--train <dir>`: directory with tensor files.
    - `--vocab <file>`: vocab file.
    - `--save_dir <dir>`: where to write checkpoints.
    - Common hyperparams: `--batch_size`, `--epoch`, `--warmup`, `--step_beta`, `--kl_anneal_iter`, `--anneal_iter`, `--print_iter`, `--save_iter`.

- `generate.py`
    - `--vocab <file>`: vocab file.
    - `--model <ckpt>`: checkpoint path.
    - `--nsamples <int>`: number of molecules to generate.

- `custom_scripts/dataset_cleaning/canonicalize_smiles.py`
    - `--input/-i <file>`: SMILES input.
    - `--output/-o <file>`: canonicalized, deduped output.

- `custom_scripts/computational_screening/tanimoto.py`
    - `--input/-i <file>`: generated SMILES (one per line).
    - `--training/-t <file>`: training/reference SMILES.
    - `--output/-o <dir>`: output directory (created if missing).
    - `--threshold <float>`: novelty cutoff (default 0.7). Outputs all CSV, threshold CSV, and novel SMILES txt.

- `custom_scripts/computational_screening/ro5.py`
    - `--input/-i <csv>`: CSV from Tanimoto (needs canonical_smiles column).
    - `--output/-o <dir>`: output directory. Writes `lipinski_all.csv`, `lipinski_drug_like.csv`, `drug_like_canonical_smiles.txt`.

- `custom_scripts/DeepPurpose/screening.py`
    - `--smiles/-s <txt>`: SMILES list (e.g., drug_like_canonical_smiles.txt).
    - `--targets/-t <csv|txt>`: protein targets (CSV with Sequence and name column, or txt with `Name: SEQ`).
    - `--output/-o <prefix>`: output prefix under results/dti_screening/ (by default).
    - `--model <name>`: DeepPurpose pretrained model (default MPNN_CNN_BindingDB).

# Visualizations Guide

This folder contains utilities to visualize screening outputs, novelty, Ro5 compliance, and molecule images.

## Prerequisites
- Python 3.8+
- Install deps (from repo root):
  ```bash
  pip install pandas matplotlib seaborn rdkit-pypi numpy
  ```

## Run Location
Change into the scripts folder so relative paths resolve:
```bash
cd visualizations/scripts
```

## Scripts
- `props.py` — Compute Ro5 props from a screening CSV and write CSV summary.
  - Input columns: `Rank`, `SMILES`.
  - Example:
    ```bash
    python props.py -i ../../results/dti_screening/screening_multi_target_results.csv -o ../../results/props/props_screening_multi_target_results.csv
    ```

- `novelty.py` — Histogram of Tanimoto similarity scores with a novelty threshold line.
  - Input column: `tanimoto_score`.
  - Example:
    ```bash
    python novelty.py -i ../results/graphs/tanimoto/tanimoto_all.csv -t 0.07 -o ../results/graphs/tanimoto/novelty_hist.png
    ```

- `ro5_bubble_graph.py` — Bubble plot encoding MW (x), LogP (y), HBD (size), HBA (color), violations (edge color).
  - Input column: `canonical_smiles` (override with `--smiles-column` if different).
  - Example:
    ```bash
    python ro5_bubble_graph.py -i ../../results/dti_screening/screening_multi_target_results.csv -o ../results/graphs/ro5/ro5_bubble.png
    ```

- `smile_visuals.py` — Generates molecule PNGs, ranking bar plot, and consistency violin plot for multi-target screening.
  - Required columns: `Rank`, `SMILES`, `Mean_Affinity`, `Std_Affinity`, and per-target columns prefixed `A0A`.
  - Example:
    ```bash
    python smile_visuals.py -i ../../results/dti_screening/screening_multi_target_results.csv \
      --output-dir ../results/molecular_images/DeepPurpose \
      --ranking-plot ../results/graphs/DeepPurpose/rank.png \
      --consistency-plot ../results/graphs/DeepPurpose/consistency.png
    ```

- `DeepPurpose/visulize_results.py` — Visualization helper for DeepPurpose outputs (ranking + consistency + molecule images).
  - Example:
    ```bash
    python DeepPurpose/visulize_results.py \
      --input ../../results/dti_screening/screening_multi_target_results.csv \
      --output ../results/molecular_images/DeepPurpose/ \
      --ranking-plot ../results/graphs/DeepPurpose/rank.png \
      --consistency-plot ../results/graphs/DeepPurpose/consistency.png
    ```

## Tips
- Paths above assume you run from `visualizations/scripts`. Adjust `..` prefixes if you run elsewhere.
- Outputs are written under `visualizations/results/` by default; feel free to override `-o/--output*` flags.
- For large datasets, prefer saving figures (`-o path.png`) instead of displaying interactively.

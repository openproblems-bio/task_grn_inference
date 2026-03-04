# GRN Evaluation in GeneRNBI

---

## Metadata

**Short Description**: How GRN predictions are evaluated in GeneRNBI — output format requirements, evaluation metrics, and how to run the evaluation pipeline.

**Authors**: GeneRNBI Team

**Version**: 1.0

**Last Updated**: March 2026

**License**: CC BY 4.0

**Commercial Use**: ✅ Allowed

---

## Overview

GeneRNBI evaluates the top 50,000 TF-gene edges (ranked by weight) from your inferred GRN against multiple reference datasets and metrics. Your output must be an AnnData `.h5ad` file in a specific format.

---

## Required Output Format

Your inferred GRN must be saved as an AnnData object with a `prediction` DataFrame in `uns`:

```python
import anndata as ad
import pandas as pd

# net must have these three columns exactly:
net = pd.DataFrame({
    'source': ['TF1', 'TF2', ...],   # transcription factor name
    'target': ['GeneA', 'GeneB', ...],  # target gene name
    'weight': [0.95, 0.87, ...]        # regulatory score (higher = stronger)
})

# weight must be cast to string
net['weight'] = net['weight'].astype(str)

output = ad.AnnData(
    X=None,
    uns={
        'method_id': 'your_method_name',   # e.g. 'pearson_corr'
        'dataset_id': 'op',                # dataset used for inference
        'prediction': net[['source', 'target', 'weight']]
    }
)
output.write('output/prediction.h5ad')
```

**Rules:**
- Only TF-gene pairs are evaluated — always filter `source` to TFs in `tf_all`
- Top 50,000 edges by absolute weight are used — sort by `|weight|` descending
- `weight` must be stored as string in the AnnData

---

## Running Evaluation (without Docker)

```bash
bash src/metrics/all_metrics/run_local.sh \
  --dataset op \
  --prediction output/prediction.h5ad \
  --score output/score.h5ad \
  --num_workers 4
```

Supported datasets: `op`, `norman`, `replogle`, and others (see `resources/grn_benchmark/`).

---

## Running Evaluation (with Docker / full pipeline)

```bash
viash run src/metrics/all_metrics/config.vsh.yaml -- \
  --dataset op \
  --prediction output/prediction.h5ad \
  --score output/score.h5ad
```

---

## Available Datasets

| Dataset | Type | Notes |
|---------|------|-------|
| `op` | Perturbation (pseudobulk) | Main benchmark dataset |
| `norman` | Perturbation | CRISPR screen |
| `replogle` | Perturbation | Large-scale screen |

Download with:
```bash
aws s3 sync s3://openproblems-data/resources/grn/grn_benchmark resources/grn_benchmark/ --no-sign-request
# Test data only:
aws s3 sync s3://openproblems-data/resources_test/grn/grn_benchmark resources_test/grn_benchmark --no-sign-request
```

---

## Evaluation Metrics

Metrics are grouped into categories. Not all metrics apply to all datasets.

- **AUROC / AUPRC**: Rank-based metrics against curated TF-target reference networks (e.g. CollecTRI)
- **Regression-based**: How well TF weights predict target gene expression changes under perturbation
- **Stability**: Consistency across subsampled datasets
- **Context specificity**: Whether the network captures cell-type-specific regulation

See `src/metrics/` for individual metric implementations.

---

## Checking Your Score Output

```python
import anndata as ad
score = ad.read_h5ad('output/score.h5ad')
print(score.uns)  # contains all metric scores
```

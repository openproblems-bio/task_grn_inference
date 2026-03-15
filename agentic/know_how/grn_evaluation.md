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
  --num_workers 4 \
  --layer lognorm \
  --reg_type ridge
```

Supported `--dataset` values: see table below. `--layer` (default: `lognorm`) and `--reg_type` (default: `ridge`) are optional.

---

## Running Evaluation (full pipeline)

```bash
bash scripts/run_grn_evaluation.sh \
  --prediction=output/prediction.h5ad \
  --save_dir=output/ \
  --dataset=op \
  --build_images=false \
  --run_local=true
```

---

## Available Datasets

All dataset names accepted by `--dataset`:

| Dataset | Cell type |
|---------|-----------|
| `op` | PBMC |
| `parsebioscience` | PBMC |
| `300BCG` | PBMC |
| `ibd_uc` | PBMC |
| `ibd_cd` | PBMC |
| `replogle` | K562 |
| `norman` | K562 |
| `xaira_HEK293T` | HEK293T |
| `xaira_HCT116` | HCT116 |
| `nakatake` | — |

Download with:
```bash
aws s3 sync s3://openproblems-data/resources/grn/grn_benchmark resources/grn_benchmark/ --no-sign-request
# Test data only:
aws s3 sync s3://openproblems-data/resources_test/grn/grn_benchmark resources_test/grn_benchmark --no-sign-request
```

---

## Evaluation Metrics

Metric names used in `config.env` and score output:

| Metric key | Description |
|------------|-------------|
| `regression` | How well TF weights predict perturbation response |
| `gs_recovery` | Recovery of gold-standard TF-target pairs |
| `tf_binding` | Overlap with ChIP-seq / UniBind / ChIP-Atlas / ReMap |
| `ws_distance` | Wasserstein distance-based score |
| `vc` | Variance consistency across datasets |
| `sem` | Structural equation model score |
| `replicate_consistency` | Consistency across biological replicates |
| `tf_recovery` | Recovery of known TF regulatory programs |

Final score metrics: `r_precision`, `r_recall`, `vc`, `sem`, `ws_precision`, `ws_recall`.
Not all metrics apply to all datasets — see `src/utils/config.env` for per-dataset metric lists.

See `src/metrics/` for individual metric implementations.

---

## Checking Your Score Output

```python
import anndata as ad
score = ad.read_h5ad('output/score.h5ad')
print(score.uns)  # contains all metric scores
```

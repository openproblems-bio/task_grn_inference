# GRN Output Format and Evaluation

---

## Metadata

**Short Description**: Specification for the GRN prediction output format (AnnData), evaluation metrics, and how to run evaluation in GeneRNBI.

**Authors**: GeneRNBI Team

**Version**: 1.0

**Last Updated**: March 2026

**License**: CC BY 4.0

**Commercial Use**: ✅ Allowed

---

## GRN Prediction Output Format

All GRN methods must output an AnnData `.h5ad` file with the following structure:

```python
import anndata as ad
import pandas as pd

# net must be a DataFrame with exactly these three columns:
net = pd.DataFrame({
    'source': ['TF1', 'TF2', ...],   # Transcription factor gene name
    'target': ['GENE1', 'GENE2', ...], # Target gene name
    'weight': [0.95, 0.87, ...]       # Regulatory score (higher = stronger)
})

# IMPORTANT: weight must be stored as string
net['weight'] = net['weight'].astype(str)

output = ad.AnnData(
    X=None,
    uns={
        'method_id': 'your_method_name',           # e.g. 'pearson_corr'
        'dataset_id': rna.uns.get('dataset_id', 'unknown'),  # e.g. 'op'
        'prediction': net[['source', 'target', 'weight']]
    }
)
output.write(par['prediction'])
```

### Key rules:
- Only **top 50,000 TF-gene edges** are evaluated — always filter/sort by weight descending
- Only **TF-gene pairs** are scored — filter `source` to known TFs using `tf_all.csv`
- `source` must be a TF from the `tf_all.csv` list
- `target` can be any gene in the expression data
- No self-loops (`source != target`)

---

## Evaluation Metrics

GeneRNBI uses multiple metrics organized in three categories. Not all metrics apply to all datasets.

| Category | Metric | Description |
|----------|--------|-------------|
| **Regression** | `regression_1` | How well TF activity predicts perturbation response |
| **Regression** | `regression_2` | Variant of regression with different normalization |
| **Classification** | `auroc` | Area under ROC for ChIP-seq validation |
| **Classification** | `auprc` | Area under precision-recall curve |
| **Network** | `connectedness` | Structural network properties |

---

## Running Evaluation

### Without Docker (recommended for testing):
```bash
bash src/metrics/all_metrics/run_local.sh \
    --dataset op \
    --prediction output/net.h5ad \
    --score output/score.h5ad \
    --num_workers 4
```

### With Docker (full pipeline):
```bash
viash run src/metrics/all_metrics/config.vsh.yaml -- \
    --dataset op \
    --prediction output/net.h5ad \
    --score output/score.h5ad
```

### Available datasets for evaluation:
- `op` — Open Problems perturbation dataset
- `norman` — Norman et al. CRISPR screen
- `replogle` — Replogle et al. genome-wide Perturb-seq

---

## Reading Evaluation Results

```python
import anndata as ad
score = ad.read_h5ad('output/score.h5ad')
print(score.uns)  # contains all metric scores as a dict
```

---

## Example: Checking Your Output is Valid

```python
import anndata as ad
import numpy as np

pred = ad.read_h5ad('output/prediction.h5ad')

# Check required fields
assert 'method_id' in pred.uns
assert 'dataset_id' in pred.uns
assert 'prediction' in pred.uns

net = pred.uns['prediction']
assert set(net.columns) >= {'source', 'target', 'weight'}
assert len(net) <= 50000, f"Too many edges: {len(net)}"
assert net['weight'].dtype == object, "weight must be string dtype"

print(f"Valid prediction: {len(net)} edges")
print(net.head())
```

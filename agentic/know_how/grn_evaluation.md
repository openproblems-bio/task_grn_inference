# GRN Evaluation in geneRNBI

---

How GRN predictions are evaluated in geneRNBI — output format requirements, evaluation metrics, and how to run the evaluation pipeline.

## Overview

Metrics:

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

Not all metrics apply to all datasets.

| Dataset | Metrics |
|---------|---------|
| `op` | regression, vc, replicate_consistency, tf_binding, sem, gs_recovery |
| `parsebioscience` | regression, vc, replicate_consistency, tf_binding, sem, gs_recovery |
| `300BCG` | regression, vc, replicate_consistency, tf_binding, sem, gs_recovery |
| `ibd_uc` | regression, tf_binding, gs_recovery, replicate_consistency |
| `ibd_cd` | regression, tf_binding, gs_recovery, replicate_consistency |
| `replogle` | regression, ws_distance, tf_recovery, tf_binding, sem, gs_recovery, vc |
| `norman` | regression, ws_distance, tf_binding, gs_recovery, vc |
| `nakatake` | regression, sem, gs_recovery, vc |
| `xaira_HEK293T` | regression, ws_distance, tf_recovery, tf_binding, sem, gs_recovery, vc |
| `xaira_HCT116` | regression, ws_distance, tf_recovery, tf_binding, sem, gs_recovery, vc |


See `src/metrics/` for individual metric implementations.


geneRNBI evaluates the top 50,000 TF-gene edges (ranked by weight) from your inferred GRN against multiple reference datasets and metrics. Your output must be an AnnData `.h5ad` file in a specific format.

---

## Required GRN model format

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


### no sub cell type within inferred GRN
This is applicable only for inferred GRNs that are not cell type specific -> for example, one GRN for PBMC.

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
### sub cell type within inferred GRN
This is applicable only for inferred GRNs that are cell type specific -> for example, one GRN for each subtype of PBMC, e.g. CD8T, CD4T.

```python
import anndata as ad
import numpy as np

pred = ad.read_h5ad('output/prediction.h5ad')

# Check required fields
assert 'method_id' in pred.uns
assert 'dataset_id' in pred.uns
assert 'prediction' in pred.uns

net = pred.uns['prediction']
assert set(net.columns) >= {'source', 'target', 'weight', 'cell_type'}
assert net.groupby(['source', 'target', 'weight']...###TODO: fix this > 50k)
assert net['weight'].dtype == object, "weight must be string dtype"

print(f"Valid prediction: {len(net)} edges")
print(net.head())


```

---

## Running Evaluation (without Docker)

```bash
bash src/metrics/all_metrics/run_local.sh \
  --dataset {dataset} \
  --prediction {grn_model}.h5ad \
  --score output/score.h5ad \
  --num_workers 4 
```

---

## Running Evaluation (full pipeline)
This needs Viash and Docker installed.

```bash
bash scripts/run_grn_evaluation.sh \
  --prediction=output/prediction.h5ad \
  --save_dir=output
  --dataset={dataset} 
  --build_images=false  # needs to be true at first run
  --run_local=true # whether to run on aws or locally
```

## Update the consensus for a given new model
If you have a new GRN model, run the following command to update the consensus files before running the evaluations. This will generally works in favour of the model as it includes GRN topology of the given GRN in the evaluation.

```bash
bash scripts/prior/run_consensus.sh --dataset {dataset} --new_model {new_grn_model_file.h5ad}
```


## Checking Your Score Output

```python
import anndata as ad
score = ad.read_h5ad('output/score.h5ad')
print(score.uns)  # contains all metric scores
```

## Comparing your model's performance to previous results
The previous benchmark scores are in resources/results/all_scores.csv for all methods, datasets, and metrics.
Of note, many metrics are there but not actually used as the final evaluation metris (we use those in stability analysis). Use src/utils/config.py -> METRICS to take out relevant metrics.

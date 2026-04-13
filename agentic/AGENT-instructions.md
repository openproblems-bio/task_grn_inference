# geneRNBI Agent Instructions

Agent-facing guide for working in the geneRNBI GRN benchmark repository.

---

## What this repo does

**geneRNBI** benchmarks GRN (gene regulatory network) inference methods across multiple datasets and evaluation metrics. The pipeline:
1. Runs GRN inference methods on OMICs data (`src/methods/`)
2. Evaluates the resulting networks against reference datasets (`src/metrics/`)
3. Scores and ranks methods on a leaderboard

---

## Key directory structure

```
src/
  methods/          # GRN inference methods (one folder per method)
  metrics/          # Evaluation metrics (one folder per metric)
  api/              # Shared Viash schemas (comp_method.yaml, comp_control_method.yaml, ...)
  utils/util.py     # Shared Python helpers (corr_net, manage_layer, parse_args, ...)

resources/
  grn_benchmark/
    inference_data/   # Full-size RNA h5ad files for inference
    evaluation_data/  # Bulk/pseudobulk/DE datasets for scoring
    prior/            # tf_all.csv, consensus priors, skeleton

resources_test/
  grn_benchmark/      # Small copies of the above — use for local development & viash test

scripts/
  run_grn_inference.sh     # Run a method on a dataset using Nextflow/Docker/Viash
  run_grn_evaluation.sh    # Run all metrics on a prediction file Nextflow/Docker/Viash
  prior/run_consensus.sh   # Rebuild consensus prior (needed for new GRN models)

scripts/local_workflows
  run_grn_inference.sh     # Run a method on a dataset using conda env
  run_grn_evaluation.sh    # Run all metrics on a prediction fileusing conda env

agentic/know_how/          # Detailed how-to guides (see below)
```

---

## Know-how guides

Read the relevant guide before starting a task:

| File | When to read it |
|------|----------------|
| `know_how/add_grn_method.md` | Adding a new GRN inference method |
| `know_how/grn_evaluation.md` | Running evaluation, understanding metrics, checking output format, comparing new scores to available scores |
| `know_how/datasets.md` | Downloading datasets, understanding dataset IDs and modalities |

---

## How to calculate overall performance of GRN models
A two-step normalization is applied to Raw scores (resources/results/all_scores.csv) to get normalized scores: first, each metric is min-max normalized to [0, 1] separately within each dataset (negative values clipped to 0); then the per-dataset scores are averaged across datasets and re-normalized to [0, 1] across models. An analogous aggregation is computed per dataset (averaging across metrics). The overall_score is then the median of these doubly-normalized values, restricted to metrics that passed applicability criteria (resources/results/metrics_kept_per_dataset.yaml). The final ranked summary table is saved to resources/results/summary.tsv.

function: src/utils/util -> compute_overall_scores

## Designing new GRN inference method
While the overall performance is important, consider also the computational costs given at resources/results/summary.tsv.
Also, remove positive control model as its good performance comes from using both inference and evaluation datasets. 

## What is consensus used in the evaluations?
To standardize the benchmark with respect to the topology discrepency between GRN models, define consensus to remove the topology confounders, such as GRN density, which could affect metrics such as Regression.

## GRN inference methods
We have two categories of GRN inference methods -> transcriptomics-based (S) and multiomics (M) based. Three datasets support M while 7 supports S. The M methods that are not applicable to  transcriptomics alone datasets are simply excluded from the benchmark (without penalizing).


## Config file
`src/utils/config.py` -> lists available datasets, datasets to metrics applicability, integrated methods, metrics used for evaluations, etc.

## Shared utilities (`src/utils/util.py`)


Import in scripts as:
```python
try:
    sys.path.append(meta['resources_dir'])
except:
    sys.path.append('src/utils')
from util import corr_net, parse_args, manage_layer
```

Key functions:

| Function | Purpose |
|----------|---------|
| `manage_layer(adata, par)` | Selects the right expression layer (`lognorm`, `X_norm`) |
| `parse_args(par)` | Merges CLI overrides into the `par` dict |
| `corr_net(adata, tf_all, par)` | Computes correlation-based GRN, filters to TFs, trims to `max_n_links` |
| `format_save_score(results_df, method_id, dataset_id, score_file)` | Saves metric scores in the standard h5ad format |

---


## Common errors

| Error | Cause | Fix |
|-------|-------|-----|
| `__merge__ file not found` | Wrong working directory | Run viash from repo root |
| `Cannot connect to Docker daemon` | Docker not running | `sudo systemctl start docker` or start Docker Desktop |
| `KeyError: dataset_id` | Missing metadata in h5ad | Use `rna.uns.get('dataset_id', 'unknown')` |
| `Empty prediction after TF filter` | Gene names don't match tf_all | Print `np.intersect1d(tf_all, rna.var_names)` to debug |
| `write() missing argument` | Used `output.write()` | Replace with `output.write_h5ad()` |
| `SyntaxError` in generated script | Used Python `open()` to write `.py` file | Use bash heredoc instead |
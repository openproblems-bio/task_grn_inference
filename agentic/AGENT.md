# geneRNBI Agent Instructions

Agent-facing guide for the geneRNBI GRN benchmark. Read this file first, then follow the task-routing table to the relevant know-how guide.

---

## 1. Framework in 60 seconds

**geneRNBI** benchmarks gene regulatory network (GRN) inference methods. The pipeline has three stages:

```
resources/grn_benchmark/inference_data/{dataset}_rna.h5ad
          [+ optional ATAC for multiomics methods]
                        │
                        ▼ (1) INFERENCE
          src/methods/{method}/script.py
                        │
                        ▼
          resources/results/{dataset}/{dataset}.{method}.prediction.h5ad
                        │
                        ▼ (2) EVALUATION
          src/metrics/all_metrics/  (or scripts/run_grn_evaluation.sh)
                        │
                        ▼
          resources/results/{dataset}/scores/{dataset}__{method}_score.h5ad
                        │
                        ▼ (3) POST-PROCESSING  (scripts/benchmark/)
          aggregate_local_scores.py   →  resources/results/benchmark/all_scores.csv
          evaluate_metric_applicability.py   →  valid metrics per dataset
          create_overview_figure.py   →  leaderboard overview plot
          plot_raw_scores.py   →  per-dataset score heatmaps
```

- **Inference**: each method takes RNA (+ optional ATAC) and outputs a ranked TF→gene edge list.
- **Evaluation**: metrics score the prediction against perturbation/ChIP-seq reference data.
- **Post-processing**: aggregates per-method score files, filters metrics by statistical validity, and generates figures. Run after all evaluations are complete.
- Methods are either **transcriptomics (S)** or **multiomics (M)**. M methods that lack ATAC for a dataset are silently skipped — not penalized.

---

## 2. Decide your execution environment first

Every command differs depending on which environment you are in. Determine this before writing or running anything.

| Environment | When to use | Entry points |
|-------------|-------------|--------------|
| **Docker / Viash** | CI, reproducible runs, adding new methods officially | `viash run src/methods/{m}/config.vsh.yaml`, `scripts/run_grn_inference.sh`, `scripts/run_grn_evaluation.sh` |
| **Conda (`genernbi`)** | Local dev, lightweight methods, all evaluation, all post-processing | `src/methods/{m}/run_local.sh`, `src/metrics/all_metrics/run_local.sh`, `scripts/local_workflows/`, `scripts/benchmark/` |
| **Singularity** | Local, heavy methods (GRNBoost2, scPRINT, Scenic/Scenic+, CellOracle, etc.) | `src/methods/{m}/run_local.sh` (these scripts internally invoke singularity) |
| **SLURM (`sbatch`)** | Cluster runs | `scripts/local_workflows/run_grn_inference.sh` (sets `run_prefix=sbatch`) |

**Local development rule:** always use `resources_test/` paths (380 MB) to verify your changes before running on full data (`resources/`).

---

## 3. Task routing

| You need to… | Read |
|--------------|------|
| Add a new GRN inference method | `agentic/know_how/add_grn_method.md` |
| Run GRN inference for a method/dataset | `agentic/know_how/add_grn_method.md` → Step 5 |
| Run evaluation or understand metric output | `agentic/know_how/grn_evaluation.md` |
| Aggregate score files into all_scores.csv | `agentic/know_how/grn_evaluation.md` → "Post-evaluation analysis" |
| Compare a new model to existing benchmark scores | `agentic/know_how/grn_evaluation.md` → "Comparing your model" |
| Regenerate leaderboard figures | `agentic/know_how/grn_evaluation.md` → "Post-evaluation analysis" |
| Understand datasets, download data | `agentic/know_how/datasets.md` |
| Check which metrics apply to which dataset | `src/utils/config.py` → `DATASETS_METRICS` |
| Check which methods are integrated | `src/utils/config.py` → `METHODS` |
| Debug a failing run | Section 7 of this file |

---

## 4. Data contracts

### 4a. Inference input

| File | Path pattern | Required fields |
|------|-------------|-----------------|
| RNA expression | `resources/grn_benchmark/inference_data/{dataset}_rna.h5ad` | `layers['lognorm']` or `layers['X_norm']`, `uns['dataset_id']`, `var_names` = gene names |
| ATAC (multiomics only) | `resources/grn_benchmark/inference_data/{dataset}_atac.h5ad` | `uns['dataset_id']` |
| TF list | `resources/grn_benchmark/prior/tf_all.csv` | one TF name per line |

### 4b. Prediction format (output of inference, input to evaluation)

```python
# uns['prediction'] must be a DataFrame with exactly these columns:
net = pd.DataFrame({'source': [...],   # TF name — must be in tf_all
                    'target': [...],   # target gene name
                    'weight': [...]})  # regulatory score, higher = stronger
net['weight'] = net['weight'].astype(str)   # weight stored as string

output = ad.AnnData(X=None, uns={
    'method_id': 'your_method',
    'dataset_id': rna.uns.get('dataset_id', 'unknown'),
    'prediction': net[['source', 'target', 'weight']]
})
output.write_h5ad(par['prediction'])
```

**Hard rules:**
- Only TF→gene pairs are evaluated — filter `source` to TFs in `tf_all` before saving.
- Top 50,000 edges by `|weight|` are used. Sort and trim before saving.
- Use `.write_h5ad()`, not `.write()`.
- Cell-type-specific GRNs: add a `cell_type` column to `net`.

### 4c. Score format (output of evaluation)

```python
score = ad.read_h5ad('output/score.h5ad')
# score.uns contains: dataset_id, method_id, metric_ids, metric_values
```

Score files are written to `resources/results/{dataset}/scores/{dataset}__{method}_score.h5ad`.

### 4d. Aggregated results (output of post-processing)

`resources/results/benchmark/all_scores.csv` — wide-format table with one row per (dataset, method) and one column per metric. This is the source for leaderboard rankings. Produced by `scripts/benchmark/aggregate_local_scores.py`.

---

## 5. Config as source of truth

**`src/utils/config.py`** — consult this before hardcoding any dataset/metric/method names.

| Variable | What it contains |
|----------|-----------------|
| `DATASETS` | Active benchmark dataset IDs |
| `DATASET_INFO` | Cell type, perturbation type, modality per dataset |
| `DATASETS_METRICS` | Which metrics apply to each dataset |
| `METRICS` | Final metric list used for scoring/ranking |
| `METRIC_THRESHOLDS` | Minimum meaningful value per metric (used by `evaluate_metric_applicability.py`) |
| `METHODS` | Integrated method IDs |
| `surrogate_names` | Display names for methods, metrics, datasets |

**`src/utils/util.py`** — shared helpers. Import pattern:

```python
try:
    sys.path.append(meta['resources_dir'])
except:
    sys.path.append('src/utils')
from util import corr_net, parse_args, manage_layer
```

| Function | Purpose |
|----------|---------|
| `manage_layer(adata, par)` | Selects the right expression layer (`lognorm` / `X_norm`) |
| `parse_args(par)` | Merges CLI overrides into the `par` dict |
| `corr_net(adata, tf_all, par)` | Correlation-based GRN, TF-filtered, trimmed to `max_n_links` |
| `format_save_score(results_df, method_id, dataset_id, score_file)` | Saves metric scores in standard h5ad format |
| `compute_overall_scores(df)` | Aggregates per-metric scores into an overall ranking |

---

## 6. Consensus prior

The **consensus prior** standardizes evaluation across GRN methods by removing topology confounders (e.g., edge density differences that would otherwise distort regression scores).

**When to update it:** only when evaluating a new GRN model not already in geneRNBI. Run this before evaluation:

```bash
bash scripts/prior/run_consensus.sh --dataset {dataset} --new_model output/net.h5ad
```

If you skip this step, the new model is evaluated without its topology being represented in the consensus, which disadvantages it on regression-based metrics.

---

## 7. Troubleshooting

Work through these steps in order before consulting the error table.

**Step 1 — Verify environment**
- Are you in the right conda env (`conda activate genernbi`) or using the correct singularity image?
- Is Docker running if using Viash? (`docker info`)
- Are you running from the **repo root**? Viash and relative paths require this.

**Step 2 — Verify paths**
- Do the input files exist? (`ls resources/grn_benchmark/inference_data/`)
- For test runs, are you using `resources_test/` paths in the `## VIASH START` block?
- Is the output directory writable and does it exist?

**Step 3 — Verify data format**
- Does `pred.uns['prediction']` have columns `source`, `target`, `weight`?
- Is `weight` dtype `object` (string)?
- Is `len(net) <= 50000`?
- Does `np.intersect1d(tf_all, rna.var_names)` return a non-empty array?

**Common errors**

| Error | Cause | Fix |
|-------|-------|-----|
| `__merge__ file not found` | Wrong working directory | Run viash from repo root |
| `Cannot connect to Docker daemon` | Docker not running | `sudo systemctl start docker` or start Docker Desktop |
| `KeyError: dataset_id` | Missing metadata in h5ad | Use `rna.uns.get('dataset_id', 'unknown')` |
| `Empty prediction after TF filter` | Gene names don't match tf_all | Print `np.intersect1d(tf_all, rna.var_names)` to debug |
| `write() missing argument` | Used `output.write()` | Replace with `output.write_h5ad()` |
| `SyntaxError` in generated script | Used Python `open()` to write `.py` file | Use bash heredoc: `cat > script.py << 'PYEOF'` |
| Method skipped silently on a dataset | ATAC file missing for a multiomics method | Expected — M methods require ATAC; check `resources/grn_benchmark/inference_data/{dataset}_atac.h5ad` |
| Score NaN or below threshold | Too few edges, or TF filter removed all rows | Check edge count and TF coverage; see `METRIC_THRESHOLDS` in config.py |
| `all_scores.csv` missing or stale | Post-processing not yet run | Run `scripts/benchmark/aggregate_local_scores.py` after all evaluations complete |

---

## 8. Key directory map

```
src/
  methods/{method}/         # One folder per inference method
    config.vsh.yaml         # Viash component config (Docker)
    script.py               # Inference logic
    run_local.sh            # Local/singularity runner
  metrics/{metric}/         # One folder per evaluation metric
  metrics/all_metrics/
    run_local.sh            # Runs all applicable metrics for a dataset
  api/                      # Shared Viash schemas
  utils/
    config.py               # Source of truth: datasets, methods, metrics
    util.py                 # Shared Python helpers

resources/grn_benchmark/
  inference_data/           # Full-size RNA/ATAC h5ad files
  evaluation_data/          # Perturbation/bulk/DE reference datasets
  prior/                    # tf_all.csv, consensus priors
resources/results/
  {dataset}/scores/         # Per-method score h5ad files
  all_scores.csv            # Aggregated leaderboard data

resources_test/grn_benchmark/   # Small copies — use for local dev & viash test

scripts/
  run_grn_inference.sh      # Full pipeline: inference (Docker/Nextflow/AWS)
  run_grn_evaluation.sh     # Full pipeline: evaluation (Docker/Nextflow/AWS)
  prior/run_consensus.sh    # Rebuild consensus prior
  local_workflows/
    run_grn_inference.sh    # Inference without Docker (conda/singularity/sbatch)
    run_grn_evaluation.sh   # Evaluation without Docker
  benchmark/
    aggregate_local_scores.py        # Collect score h5ad files → all_scores.csv
    evaluate_metric_applicability.py # Filter metrics by CV and threshold validity
    create_overview_figure.py        # Generate leaderboard overview plot
    plot_raw_scores.py               # Generate per-dataset score heatmaps

agentic/know_how/           # Detailed how-to guides (see task routing table above)
```

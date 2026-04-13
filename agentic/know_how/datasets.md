# geneRNBI Datasets


There are **10 integrated datasets** in the benchmark.

| Dataset | Display Name | Cell Type | Perturbation Type | Inference Data | Modality |
|---------|-------------|-----------|-------------------|----------------|----------|
| `op` | OPSCA | PBMC | Drugs | sc | Transcriptomics |
| `parsebioscience` | ParseBioscience | PBMC | Cytokines | sc/bulk | Transcriptomics |
| `300BCG` | 300BCG | PBMC | Chemicals | sc | Transcriptomics |
| `ibd_uc` | IBD:UC | PBMC | Chemicals/bacteria | sc | Multiomics |
| `ibd_cd` | IBD:CD | PBMC | Chemicals/bacteria | sc | Multiomics |
| `replogle` | Replogle | K562 | Knockout | sc/bulk | Transcriptomics |
| `norman` | Norman | K562 | Activation | sc | Transcriptomics |
| `nakatake` | Nakatake | — | — | sc/bulk | Transcriptomics |
| `xaira_HEK293T` | Xaira:HEK293T | HEK293T | Knockout | sc/bulk | Transcriptomics |
| `xaira_HCT116` | Xaira:HCT116 | HCT116 | Knockout | sc/bulk | Transcriptomics |


## Overview

All datasets are hosted on S3 and require `awscli` to download. Install it with:

```bash
pip install awscli
```

All downloads are public — no AWS credentials needed (`--no-sign-request`).

---

## Main benchmark datasets

Inference data + evaluation data + prior files (TF lists, consensus priors):

```bash
aws s3 sync s3://openproblems-data/resources/grn/grn_benchmark resources/grn_benchmark/ --no-sign-request
```

- Inference datasets: `resources/grn_benchmark/inference_data/`
- Evaluation datasets: `resources/grn_benchmark/evaluation_data/`
- Prior files (TF lists, consensus): `resources/grn_benchmark/prior/`

## Test datasets (small, for local development)

```bash
aws s3 sync s3://openproblems-data/resources_test/grn/grn_benchmark resources_test/grn_benchmark --no-sign-request
```

Use these paths in the `## VIASH START` block during local testing:

```python
par = {
    'rna': 'resources_test/grn_benchmark/inference_data/op_rna.h5ad',
    'tf_all': 'resources_test/grn_benchmark/prior/tf_all.csv',
    'prediction': 'output/prediction.h5ad'
}
```

## Extended datasets

Large single-cell perturbation datasets (Replogle, Xaira, Parse) and full pseudobulked versions:

```bash
aws s3 sync s3://openproblems-data/resources/grn/extended_data/ resources/extended_data/ --no-sign-request
```

## Raw / unprocessed data

```bash
aws s3 sync s3://openproblems-data/resources/grn/datasets_raw/ resources/datasets_raw/ --no-sign-request
```

## Pre-computed results (leaderboard)

```bash
aws s3 sync s3://openproblems-data/resources/grn/results resources/results/ --no-sign-request
```

---

## Running evaluation without Docker

If Docker is unavailable, run all metrics locally using:

```bash
bash src/metrics/all_metrics/run_local.sh \
  --dataset <dataset_name> \
  --prediction=<inferred_grn.h5ad> \
  --score <output_score.h5ad> \
  --num_workers <N>
```

Example:

```bash
bash src/metrics/all_metrics/run_local.sh \
  --dataset op \
  --prediction=resources/grn_models/op/collectri.h5ad \
  --score=output/score.h5ad \
  --num_workers=20
```

If evaluating a new GRN model not already in geneRNBI, first generate its consensus prior:

```bash
bash scripts/prior/run_consensus.sh --dataset op --new_model {new_grn_model_file.h5ad}
```

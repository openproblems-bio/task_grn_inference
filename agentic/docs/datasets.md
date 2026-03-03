# GRN Benchmark Datasets

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

## Metrics per Dataset

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

# GeneRNBI General Reference

---

## Metadata

**Short Description**: Authoritative reference for all datasets, methods, and metrics in GeneRNBI — names, groupings, applicability, and final evaluation selection.

**Authors**: GeneRNBI Team

**Version**: 1.0

**Last Updated**: March 2026

---

## Datasets (10 active)

| ID | Display Name | Cell Type | Modality | Perturbation Type |
|----|--------------|-----------|----------|-------------------|
| op | OPSCA | PBMC | Multiomics | Drugs |
| ibd_uc | IBD:UC | PBMC | Multiomics | Chemicals/bacteria |
| ibd_cd | IBD:CD | PBMC | Multiomics | Chemicals/bacteria |
| 300BCG | 300BCG | PBMC | Transcriptomics | Chemicals |
| parsebioscience | ParseBioscience | PBMC | Transcriptomics | Cytokines |
| replogle | Replogle | K562 | Transcriptomics | Knockout |
| norman | Norman | K562 | Transcriptomics | Activation |
| nakatake | Nakatake | SEES3 (PSC) | Transcriptomics | Overexpression |
| xaira_HEK293T | Xaira:HEK293T | HEK293T | Transcriptomics | Knockout |
| xaira_HCT116 | Xaira:HCT116 | HCT116 | Transcriptomics | Knockout |

(adamson is defined in the codebase but currently excluded from active evaluation)

---

## Methods (15)

| ID | Display Name | Modality |
|----|--------------|----------|
| positive_control | Positive Ctrl | — |
| negative_control | Negative Ctrl | — |
| pearson_corr | Pearson Corr. | Transcriptomics |
| spearman_corr | Spearman Corr. | Transcriptomics |
| ppcor | PPCOR | Transcriptomics |
| grnboost | GRNBoost2 | Transcriptomics |
| portia | Portia | Transcriptomics |
| scenic | Scenic | Transcriptomics |
| scgpt | scGPT | Transcriptomics |
| geneformer | GeneFormer | Transcriptomics |
| scprint | scPRINT | Transcriptomics |
| scenicplus | Scenic+ | Multiomics |
| celloracle | CellOracle | Multiomics |
| figr | FigR | Multiomics |
| granie | GRaNIE | Multiomics |
| scglue | scGLUE | Multiomics |

---

## Metric Groups and Their Datasets

Each dataset supports a subset of metric groups:

| Metric Group | Datasets |
|---|---|
| regression | all 10 datasets |
| ws_distance | replogle, norman, xaira_HEK293T, xaira_HCT116 |
| tf_recovery | replogle, xaira_HEK293T, xaira_HCT116 |
| tf_binding | replogle, norman, op, 300BCG, ibd_uc, ibd_cd, parsebioscience, xaira_HEK293T, xaira_HCT116 |
| sem | replogle, nakatake, op, 300BCG, parsebioscience, xaira_HEK293T, xaira_HCT116 |
| gs_recovery | replogle, norman, nakatake, op, 300BCG, ibd_uc, ibd_cd, parsebioscience, xaira_HEK293T, xaira_HCT116 |
| vc | replogle, norman, nakatake, op, 300BCG, parsebioscience, xaira_HEK293T, xaira_HCT116 |
| replicate_consistency | op, 300BCG, ibd_uc, ibd_cd, parsebioscience |

---

## Fine-Grained Metrics (11 total — all computed and reported)

All 11 metrics below are computed and reported for every applicable method/dataset combination. A minimum threshold must be exceeded for a metric to be considered meaningful on a given dataset. Only a subset passed additional quality control criteria (context-specificity and robustness in stability analysis) and are used for the final ranking (see the manuscript).

| Metric ID | Display Name | Metric Group | Threshold | Final ranking |
|-----------|--------------|--------------|-----------|---------------|
| r_precision | Regression (precision) | regression | 0.1 | ✅ yes |
| r_recall | Regression (recall) | regression | 0.1 | ✅ yes |
| ws_precision | WS (precision) | ws_distance | 0.5 | ✅ yes |
| ws_recall | WS (recall) | ws_distance | 0.5 | ✅ yes |
| vc | Virtual cell | vc | 0.1 | ✅ yes |
| sem | SEM | sem | 0.1 | ✅ yes |
| t_rec_precision | TF recovery (precision) | tf_recovery | 2.0 | ❌ no |
| t_rec_recall | TF recovery (recall) | tf_recovery | 2.0 | ❌ no |
| replicate_consistency | Replicate consistency | replicate_consistency | 0.3 | ❌ no |
| tfb_f1 | TF binding | tf_binding | 0.05 | ❌ no |
| gs_f1 | Gene sets recovery | gs_recovery | 0.1 | ❌ no |

---

## Final Ranking Metrics (6 — passed quality control)

Only these 6 metrics are used in the final score / leaderboard ranking:

- r_precision
- r_recall
- vc
- sem
- ws_precision
- ws_recall

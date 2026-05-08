
Leaderboard
=================

The overall comparative performance of the integrated GRN inference methods is summarized below.
Only metrics that pass the applicability criteria for a given dataset are used to rank the methods.

The table below shows which metrics are applicable per dataset (✓ = applicable by metric requirements; ★ = also passes quality criteria).

.. image:: images/metric_applicability.png
   :width: 60%
   :align: center

----

Overall score normalization
----------------------------

The overall leaderboard score is computed in three steps:

1. **Per-dataset min-max normalization** — for each metric within each dataset, raw scores are clipped to zero for negative values and scaled to [0, 1] across methods.
2. **Percentile ranking** — methods are ranked by percentile across all applicable (dataset, metric) pairs (1 = best).
3. **Modality-balanced aggregation** — the final overall score is the mean of per-modality median percentile ranks (transcriptomics and multiomics datasets contribute equally, regardless of how many datasets exist in each modality). Methods that are not applicable to a modality are excluded from that modality's rank.

This ensures that a method with strong transcriptomics performance is not penalized simply because there are more transcriptomics datasets than multiomics ones.

:download:`Download raw benchmark scores (CSV) <data/all_scores.csv>`

----

Summary
--------

.. image:: images/metrics_map.png
   :width: 45%
   :align: center

----

.. image:: images/summary_figure.png
   :width: 100%
   :align: center

----

Raw scores per dataset
-----------------------

The heatmaps below show raw (unnormalized) metric scores for each dataset. Rows are GRN inference methods; columns are individual sub-metrics. Color represents the raw score value — higher is better for all metrics. Grey cells indicate that a metric is not applicable to that dataset or that the method did not produce a valid result.

Note that sub-metric names in these heatmaps differ from the display names used in the leaderboard. See ``surrogate_names`` in ``src/utils/config.py`` for the mapping.

Multiomics datasets
^^^^^^^^^^^^^^^^^^^^

**OPSCA** (PBMC, drug perturbations — scRNA + scATAC)

.. image:: images/raw_scores_op.png
   :width: 70%
   :align: center

----

**MSCIC** (BMMC, observational — snRNA + snATAC, 10x Multiome)

.. image:: images/raw_scores_MSCIC.png
   :width: 70%
   :align: center

----

Transcriptomics datasets
^^^^^^^^^^^^^^^^^^^^^^^^^

**Nakatake** (PSC, TF overexpression)

.. image:: images/raw_scores_nakatake.png
   :width: 70%
   :align: center

----

**Norman** (K562, CRISPRa activation)

.. image:: images/raw_scores_norman.png
   :width: 70%
   :align: center

----

**Replogle** (K562, CRISPRi knockout)

.. image:: images/raw_scores_replogle.png
   :width: 70%
   :align: center

----

**300BCG** (PBMC, chemical perturbations)

.. image:: images/raw_scores_300BCG.png
   :width: 70%
   :align: center

----

**ParseBioscience** (PBMC, cytokine stimulation)

.. image:: images/raw_scores_parsebioscience.png
   :width: 70%
   :align: center

----

**Xaira HEK293T** (HEK293T, CRISPRi knockout)

.. image:: images/raw_scores_xaira_HEK293T.png
   :width: 70%
   :align: center

----

**Xaira HCT116** (HCT116, CRISPRi knockout)

.. image:: images/raw_scores_xaira_HCT116.png
   :width: 70%
   :align: center

----

**SoundLife** (CD4 T cells, longitudinal observational)

.. image:: images/raw_scores_soundlife.png
   :width: 70%
   :align: center

----

**SoundLife: Vaccine** (B cells, flu vaccination)

.. image:: images/raw_scores_soundlife_vaccine.png
   :width: 70%
   :align: center

----

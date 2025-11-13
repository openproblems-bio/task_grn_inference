# RC_TF_ACT: Replica Consistency of TF Activity Metric

## Overview

This metric measures the consistency of transcription factor (TF) activities across biological replicates (donors). The intuition is that a good GRN should produce consistent TF activity scores across replicas within the same experimental condition (perturbation × cell type).

## Methodology

### 1. TF Activity Calculation
- Uses **Decoupler's ULM (Univariate Linear Model)** to infer TF activities from gene expression and GRN
- For each sample, calculates activity scores for all TFs based on their target genes

### 2. Grouping Strategy
- Groups samples by `perturbation × celltype × donor` (configurable via `DATASET_GROUPS[dataset_id]['rc_tf_ac']`)
- Within each `perturbation × celltype` combination, measures consistency across donors

### 3. Dispersion Metric
- For each TF and each experimental condition, calculates **Coefficient of Variation (CV)** across donors:
  - `CV = std(activities) / |mean(activities)|`
- Lower CV = more consistent TF activity = better GRN

### 4. Consistency Score
- Converts dispersion to consistency: `consistency = 1 - normalized_dispersion`
- Normalization uses robust statistics (median + 2×MAD)

### 5. Baseline Comparison
- Creates random baseline by shuffling target genes while preserving degree distribution
- Compares original GRN consistency vs random baseline

## Output Metrics

1. **rc_tf_act_precision**: Ratio of mean consistency (original / baseline)
   - Higher is better
   - Values > 1.0 indicate better than random

2. **rc_tf_act_balanced**: -log10(p-value) from Wilcoxon signed-rank test
   - Higher is better
   - Measures statistical significance of improvement over baseline

## Supported Datasets

Only datasets with `rc_tf_ac` grouping defined in `config.py`:
- `op`
- `300BCG`
- `parsebioscience`
- `ibd`

## Usage

### Local Testing
```bash
python src/metrics/rc_tf_act/script.py \
    --prediction resources/results/op/op.grnboost.grnboost.prediction.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --score output/rc_tf_act/test_score.h5ad \
    --min_targets 5
```

### Batch Processing
```bash
sbatch src/metrics/rc_tf_act/run_local.sh
```

## Parameters

- `--min_targets`: Minimum number of target genes required per TF (default: 5)
- `--layer`: Expression layer to use (default: 'lognorm')

## Dependencies

- decoupler
- scipy
- numpy
- pandas
- anndata

## Key Differences from RC Metric

- **RC (replica_consistency)**: Measures correlation consistency of target genes
- **RC_TF_ACT**: Measures TF activity score consistency across replicates
  - More biologically interpretable
  - Directly tests if GRN produces consistent regulatory signals
  - Less sensitive to edge direction permutations (uses ULM which considers edge weights)

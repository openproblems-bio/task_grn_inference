# VC Metric Robustness Analysis - Degree-Preserved Permutations

## Overview

This directory contains tools for analyzing the robustness of the VC (Variance Component) metric under different types of GRN (Gene Regulatory Network) permutations, with a focus on **degree-preserved permutations**.

## What is `create_grn_baseline`?

The `create_grn_baseline` function (in `util.py`) creates a **degree-preserved permutation** of a GRN:

- **Preserves in-degree and out-degree**: Each gene maintains the same number of regulators and targets
- **Rewires topology**: Shuffles which genes regulate which genes (changes network structure)
- **Reassigns weights**: Redistributes original edge weights to new edges deterministically
- **Removes self-loops**: Ensures no gene regulates itself

### Purpose
This is useful for testing whether the VC metric is sensitive to:
1. **Degree distribution** (how many connections each gene has)
2. **Network topology** (which specific genes are connected)

If the metric changes significantly with degree-preserved permutation, it means the metric cares about topology beyond just degree distribution.

## Files Created

### 1. `analyze_robustness_degree_preserved.py`
Performs degree-preserved permutation analysis at different intensities.

**Key features:**
- Mixes original and baseline networks at intensities: 0%, 10%, 20%, 30%, 50%, 70%, 100%
- Evaluates VC metric at each intensity level
- Creates visualization showing metric degradation
- Outputs CSV with detailed results

**Usage:**
```bash
python src/metrics/vc/analyze_robustness_degree_preserved.py \
    --prediction resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --output_dir output/vc/robustness_degree_preserved \
    --intensities 0.0 0.1 0.2 0.3 0.5 0.7 1.0
```

**Outputs:**
- `degree_preserved_robustness.csv` - Detailed metrics
- `degree_preserved_robustness_plot.png` - Visualization

### 2. `run_robustness_degree_preserved.sh`
SLURM script to run degree-preserved analysis on cluster.

**Usage:**
```bash
sbatch src/metrics/vc/run_robustness_degree_preserved.sh
```

### 3. `analyze_robustness_comparison.py`
Comprehensive script that runs ALL permutation types and compares them:
- Standard permutations (net, sign, direction)
- Degree-preserved permutations
- Creates comparative visualizations
- Generates summary report

**Usage:**
```bash
python src/metrics/vc/analyze_robustness_comparison.py \
    --prediction resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --output_dir output/vc/robustness_comprehensive
```

**Outputs:**
- `standard/robustness_analysis.csv` - Standard permutation results
- `degree_preserved/degree_preserved_robustness.csv` - Degree-preserved results
- `comparison_all_permutations.png` - Main comparison plot
- `detailed_comparison.png` - Detailed subplots
- `summary_report.txt` - Text summary with key insights

### 4. `run_robustness_comparison.sh`
SLURM script to run comprehensive comparison.

**Usage:**
```bash
sbatch src/metrics/vc/run_robustness_comparison.sh
```

## Interpretation Guide

### What Each Permutation Type Tests

| Permutation Type | What It Shuffles | What It Preserves | Tests Sensitivity To |
|-----------------|------------------|-------------------|---------------------|
| **direction** | Edge directions (TF↔target) | Weights, connections | Directionality |
| **sign** | Edge signs (+/-) | Connections, magnitudes | Activation vs repression |
| **weight** | Edge weights (adds noise) | Connections, topology | Weight precision |
| **net** | All connections randomly | Total number of edges | Network structure |
| **degree_preserved** | Topology (who→who) | In/out degrees | Topology beyond degree |

### Expected Results

If the metric is **robust to degree-preserved permutation**:
- → Metric primarily cares about **degree distribution**
- → Specific topology (which genes connect to which) matters less

If the metric **changes significantly** with degree-preserved permutation:
- → Metric is sensitive to **specific network topology**
- → The exact wiring pattern matters, not just degree

### Comparing Degree-Preserved vs. Random Net Shuffle

The key comparison is:
```
Random Net Shuffle impact  vs.  Degree-Preserved Shuffle impact
```

- If `|Δ_net| > |Δ_degree_preserved|`: Degree distribution is important
- If `|Δ_degree_preserved| ≈ |Δ_net|`: Topology matters beyond just degree

## Example Workflow

### Option 1: Run only degree-preserved analysis
```bash
# Activate environment
source ~/.bash_profile
conda activate py10

# Run analysis
python src/metrics/vc/analyze_robustness_degree_preserved.py \
    --prediction resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --output_dir output/vc/robustness_degree_preserved \
    --intensities 0.0 0.1 0.2 0.3 0.5 0.7 1.0
```

### Option 2: Run comprehensive comparison
```bash
# Activate environment
source ~/.bash_profile
conda activate py10

# Run comprehensive analysis
python src/metrics/vc/analyze_robustness_comparison.py \
    --prediction resources/results/op/op.pearson_corr.pearson_corr.prediction.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --output_dir output/vc/robustness_comprehensive
```

### Option 3: Submit to cluster
```bash
# For degree-preserved only
sbatch src/metrics/vc/run_robustness_degree_preserved.sh

# For comprehensive comparison
sbatch src/metrics/vc/run_robustness_comparison.sh
```

## Customization

### Change intensity levels
Edit the `--intensities` parameter:
```bash
--intensities 0.0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
```

### Change dataset/method
Edit the configuration in the shell scripts:
```bash
DATASET="your_dataset"
METHOD="your_method"
```

### Change evaluation settings
```bash
--n_top_genes 5000  # Use more genes
```

## Output Files Explained

### CSV Files
- `intensity`: Permutation intensity (0-100%)
- `vc`: VC metric R² score
- `n_edges_original`: Number of edges in original network
- `n_edges_permuted`: Number of edges after permutation

### Plots
- Line plots showing metric degradation
- Annotations for key points (0% and 100%)
- Shaded regions showing robustness
- Comparison plots overlaying multiple permutation types

### Summary Report
- Initial vs. final R² values
- Absolute and percent changes
- Area under curve (AUC) metric
- Key insights about topology importance

## Notes

- All analyses use the same random seed (42) for reproducibility
- Degree-preserved permutation uses deterministic mixing
- Results are automatically saved and visualized
- Temporary files are cleaned up after evaluation

## References

The `create_grn_baseline` function implements a degree-preserving network randomization algorithm, useful for null model generation and topology analysis.

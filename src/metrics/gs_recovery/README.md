# Pathway Annotation Fidelity Metric

## Overview

This metric evaluates how well a Gene Regulatory Network (GRN) represents known biology by testing whether predicted transcription factor (TF) targets are enriched for coherent biological pathways. It uses over-representation analysis (ORA) with MSigDB Hallmark pathways.

## Biological Rationale

**Key Assumption**: Good GRNs should capture biologically meaningful regulatory relationships. TFs typically regulate functionally related genes (e.g., cell cycle genes, metabolic pathways), so enrichment for known pathways suggests the GRN captures real regulatory modules.

**Why This Works**:
- TFs coordinate biological processes by regulating gene sets
- Pathway databases curate known functional gene modules
- Enrichment indicates biological coherence vs random connections

## Metric Calculation

### Algorithm

For each TF in the GRN:

1. **Extract Targets**: Get predicted target genes from the GRN
2. **Filter**: Apply target count filters (default: 10-500 targets)
3. **Pathway Enrichment**: Perform Fisher's exact test for each pathway
4. **Multiple Testing**: Apply Benjamini-Hochberg FDR correction
5. **Score Calculation**: 
   - Raw score: `-log10(best p-value)`
   - Normalized score: Account for random expectation
   - Specificity: Penalize non-specific enrichment

### Aggregation

- **pf_grn**: Mean pathway fidelity for TFs present in GRN (primary metric)
- **pf_all**: Mean for all TFs (missing TFs scored as 0)
- **pf_spec**: Mean specificity score (1 / (1 + log(n_pathways)))
- **coverage_pct**: Percentage of TFs with significant enrichment (FDR < 0.05)

## Output Metrics

| Metric | Description | Interpretation |
|--------|-------------|----------------|
| `pf_grn` | Pathway fidelity (GRN TFs) | **Primary metric**. Higher = better biological coherence |
| `pf_all` | Pathway fidelity (all TFs) | Penalizes missing TFs |
| `pf_spec` | Specificity score | Higher = more specific enrichment |
| `coverage_pct` | % TFs with enrichment | Higher = broader biological coverage |

## Parameters

### Required
- `prediction`: Path to predicted GRN (h5ad format)
- `evaluation_data`: Path to evaluation dataset (h5ad format)
- `tf_all`: Path to TF list (CSV)
- `pathway_file`: Path to GMT pathway file

### Optional
- `fdr_threshold`: FDR cutoff for significance (default: 0.05)
- `min_pathway_size`: Minimum genes per pathway (default: 5)
- `max_pathway_size`: Maximum genes per pathway (default: 500)
- `min_targets`: Minimum TF targets to evaluate (default: 10)
- `max_targets`: Maximum TF targets to use (default: 500)

## Usage

### Basic Usage

```bash
# Run the metric
viash run src/metrics/experimental/annotation/config.vsh.yaml -- \
  --prediction resources/results/op/op.method.method.prediction.h5ad \
  --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
  --tf_all resources/grn_benchmark/prior/tf_all.csv \
  --pathway_file /path/to/h.all.v2024.1.Hs.symbols.gmt \
  --score output/pathway_annotation_score.h5ad
```

### Robustness Analysis

```bash
# Run comprehensive robustness analysis
python src/metrics/experimental/annotation/analyze_robustness.py \
  --prediction resources/results/op/op.method.method.prediction.h5ad \
  --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
  --output_dir output/robustness_analysis
```

This generates:
- Parameter sensitivity analysis
- Baseline comparison with randomized GRNs
- Visualization plots

## Strengths

✅ **Biologically Interpretable**: Directly tests functional coherence
✅ **Complementary**: Different from topology-based metrics
✅ **Normalized Scoring**: Accounts for random expectation
✅ **Specificity Control**: Penalizes non-specific enrichment
✅ **Robust**: Multiple pathways tested, FDR-corrected

## Limitations & Caveats

⚠️ **Annotation Bias**: Well-studied TFs (MYC, TP53) may appear better
- *Mitigation*: Compare to randomized baselines

⚠️ **Circular Reasoning**: Pathways derived from similar data types
- *Mitigation*: Use orthogonal pathway sources when possible

⚠️ **Context Specificity**: TF-pathway relationships vary by cell type
- *Mitigation*: Use context-appropriate pathway databases

⚠️ **False Positives Can Score Well**: Random edges to hub genes
- *Mitigation*: Degree-preserving randomization baseline

⚠️ **Incomplete Annotations**: Not all genes/pathways well-annotated
- *Mitigation*: Focus on well-characterized pathways (Hallmark set)

⚠️ **Statistical Power**: Low power for TFs with few targets
- *Mitigation*: Minimum target threshold (default: 10)

## Interpretation Guide

### Score Ranges

- **pf_grn > 0.5**: Excellent biological coherence
- **pf_grn 0.3-0.5**: Good coherence
- **pf_grn 0.1-0.3**: Moderate coherence
- **pf_grn < 0.1**: Poor coherence (similar to random)

### What Makes a Good Score?

A GRN with high pathway fidelity should:
1. Have many TFs with significant enrichment (high `coverage_pct`)
2. Show specific enrichment (high `pf_spec`, not many pathways)
3. Outperform degree-preserving randomized baselines
4. Be robust to parameter changes (FDR threshold, target ranges)

## Examples

### Example 1: High-Quality GRN
```
pf_grn:          0.623
pf_all:          0.581
pf_spec:         0.742
coverage_pct:    78.3%
```
**Interpretation**: Excellent biological coherence. 78% of TFs show pathway enrichment with specific, strong signals.

### Example 2: Poor-Quality GRN
```
pf_grn:          0.087
pf_all:          0.042
pf_spec:         0.231
coverage_pct:    12.1%
```
**Interpretation**: Weak biological signal. Few TFs show enrichment, suggesting many random connections.

## Validation

The metric has been validated by:
1. **Baseline Comparison**: Real GRNs score higher than randomized networks
2. **Known Positives**: Well-studied TF-pathway pairs (e.g., MYC→cell cycle)
3. **Robustness**: Stable across parameter ranges
4. **Correlation**: Correlates with other quality metrics (precision, recall)

## References

- **MSigDB Hallmark Pathways**: Liberzon et al., Cell Syst. 2015
- **Fisher's Exact Test**: Fisher, 1922
- **FDR Correction**: Benjamini & Hochberg, JRSS 1995

## Files

- `config.vsh.yaml`: Viash configuration
- `script.py`: Main entry point
- `helper.py`: Core metric implementation
- `analyze_robustness.py`: Robustness analysis tools
- `README.md`: This documentation

## Contact

For questions or issues, please contact the GRN Benchmark team.

## Changelog

- **v1.0.0** (2025-01): Initial implementation with MSigDB Hallmark pathways

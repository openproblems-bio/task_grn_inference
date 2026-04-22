# MSCIC Processing Log

## Dataset
- **Source**: GSE194122 — NeurIPS 2021 Multimodal Single-Cell Integration Challenge
- **Tissue**: BMMC (healthy)
- **Protocol**: 10x Multiome (RNA + ATAC)
- **Donors**: 10 unique, 4 sites
- **Cells**: 69,249 total

## Split Strategy
- **Evaluation** (16,568 cells): donor 15078 @ site2 + site4, donor 18303 @ site3
  - Multi-site replicates for `replicate_consistency` metric
- **Inference** (52,681 cells): all remaining cells

## Steps
1. Copy raw file: `GSE194122_multiome_BMMC_processed.h5ad.gz` → decompressed as `resources/datasets_raw/MSCIC.h5ad`
2. Created `src/process_data/main/MSCIC/script.py`
3. Script splits RNA/ATAC modalities, applies inference/evaluation split, formats obs/var, normalises
4. Outputs:
   - `resources/grn_benchmark/inference_data/MSCIC_rna.h5ad`
   - `resources/grn_benchmark/inference_data/MSCIC_atac.h5ad`
   - `resources/grn_benchmark/evaluation_data/MSCIC_rna_all.h5ad` (pseudobulked eval RNA)

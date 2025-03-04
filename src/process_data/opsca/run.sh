viash run src/process_data/opsca/config.vsh.yaml -- \
    --op_perturbation_raw resources/datasets_raw/op_perturbation_sc_counts.h5ad \
    --op_multiome resources/datasets_raw/op_multiome_sc_counts.h5ad \
    --op_perturbation_bulk resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --op_rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --op_atac resources/grn_benchmark/inference_data/op_atac.h5ad \
    --op_rna_test resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --op_atac_test resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --op_perturbation_bulk_test resources_test/grn_benchmark/evaluation_data/op_bulk.h5ad

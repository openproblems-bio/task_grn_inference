viash run src/process_data/test_data/config.novsh.yaml -- \
    --rna resources/grn_benchmark/inference_datasets/op_rna.h5ad --rna_test resources_test/grn_benchmark/inference_datasets//op_rna.h5ad \
    --atac resources/grn_benchmark/inference_datasets/op_atac.h5ad --atac_test resources_test/grn_benchmark/inference_datasets//op_atac.h5ad \
    --perturbation_data resources/grn_benchmark/evaluation_data/op.h5ad --perturbation_data_test resources_test/grn_benchmark/evaluation_data/op.h5ad \
    --multiomics_counts resources/grn_benchmark/datasets_raw/op_multiome_sc_counts.h5ad --multiomics_counts_test resources_test/grn_benchmark/datasets_raw/op_multiome_sc_counts.h5ad \
    # --perturbation_counts resources/datasets_raw/op_perturbation_sc_counts.h5ad --perturbation_counts_test resources_test/datasets_raw/op_perturbation_sc_counts.h5ad \
    # --perturbation_counts 'resources/datasets_raw/adamson_sc_counts.h5ad' --perturbation_counts_test resources_test/datasets_raw/adamson_sc_counts.h5ad

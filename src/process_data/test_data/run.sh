viash run src/process_data/test_data/config.novsh.yaml -- \
    --rna resources/inference_datasets/op_rna.h5ad --rna_test resources_test/inference_datasets/op_rna.h5ad \
    --atac resources/inference_datasets/op_atac.h5ad --atac_test resources_test/inference_datasets/op_atac.h5ad \
    --perturbation_data resources/evaluation_datasets/op_perturbation.h5ad --perturbation_data_test resources_test/evaluation_datasets/op_perturbation.h5ad \
    --multiomics_counts resources/datasets_raw/op_multiome_sc_counts.h5ad --multiomics_counts_test resources_test/datasets_raw/op_multiome_sc_counts.h5ad \
    # --perturbation_counts resources/datasets_raw/op_perturbation_sc_counts.h5ad --perturbation_counts_test resources_test/datasets_raw/op_perturbation_sc_counts.h5ad \
    --perturbation_counts 'resources/datasets_raw/adamson_sc_counts.h5ad' --perturbation_counts_test resources_test/datasets_raw/adamson_sc_counts.h5ad

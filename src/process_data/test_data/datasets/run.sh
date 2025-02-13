viash run src/process_data/test_data/config.novsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --rna_test resources_test/grn_benchmark/inference_data//op_rna.h5ad \
    --atac resources/grn_benchmark/inference_data/op_atac.h5ad \
    --atac_test resources_test/grn_benchmark/inference_data//op_atac.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data/op.h5ad \
    --evaluation_data_test resources_test/grn_benchmark/evaluation_data/op.h5ad 
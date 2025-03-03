viash run src/methods/single_omics/scenic/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
    --prediction output/scenic_prediction.csv \
    --temp_dir output/scenic
viash run src/methods/single_omics/grnboost2/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
    --prediction output/prediction.h5ad \
    --temp_dir output/grnboost2
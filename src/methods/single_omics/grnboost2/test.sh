

viash run src/methods/single_omics/grnboost2/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --prediction output/grnboost2.h5ad
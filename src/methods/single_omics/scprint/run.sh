viash run src/methods/single_omics/scprint/config.vsh.yaml -- \
    --rna resources_test/inference_datasets/op_rna.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --prediction output/prediction.h5ad
viash run src/methods/single_omics/scgpt/config.vsh.yaml -- \
    --rna resources_test/inference_datasets/op_rna.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --prediction output/prediction.h5ad
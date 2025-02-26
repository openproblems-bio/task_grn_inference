viash run src/methods/multi_omics/granie/config.vsh.yaml -- \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --prediction output/prediction.csv \
    --tf_all resources_test/grn_benchmark/prior/tf_all.csv

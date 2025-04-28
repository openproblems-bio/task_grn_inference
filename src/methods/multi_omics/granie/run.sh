viash run src/methods/multi_omics/granie/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources/grn_benchmark/inference_data/op_atac.h5ad \
    --prediction output/prediction.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv

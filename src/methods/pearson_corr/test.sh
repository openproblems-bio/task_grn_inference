

viash run src/methods/pearson_corr/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --prediction output/pearson_net.h5ad \
    --num_workers 1
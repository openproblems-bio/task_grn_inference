viash run src/metrics/regression_1/config.vsh.yaml -- \
    --prediction output/pearson_net.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data//op_bulk.h5ad \
    --score output/score.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --num_workers 1


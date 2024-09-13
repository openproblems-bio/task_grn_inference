viash run src/metrics/regression_1/config.vsh.yaml -- \
    --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --prediction output/baseline_corr.csv \
    --score output/score.h5ad
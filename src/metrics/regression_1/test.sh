viash run src/metrics/regression_1/config.vsh.yaml -- \
    --perturbation_data resources/grn-benchmark/perturbation_data.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --prediction resources/grn_models/baselines/baseline_pearson_causal_impute.csv \
    --score output/score.h5ad
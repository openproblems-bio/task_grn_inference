viash run src/metrics/regression_1/config.vsh.yaml -- \
    --prediction resources/grn_models/norman/grnboost2.h5ad \
    --dataset_id norman --evaluation_data resources/grn_benchmark/evaluation_datasets//norman_perturbation.h5ad \
    --score output/score.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv
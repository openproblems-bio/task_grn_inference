viash run src/metrics/regression_2/config.vsh.yaml -- \
    --prediction resources/grn_models/norman/grnboost2.csv \
    --dataset_id norman --evaluation_data resources/evaluation_datasets/norman_perturbation.h5ad \
    --score output/score.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --regulators_consensus resources/prior/regulators_consensus_norman.json 
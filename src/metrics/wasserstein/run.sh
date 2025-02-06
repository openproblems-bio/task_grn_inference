viash run src/metrics/wasserstein/config.vsh.yaml -- \
    --prediction resources/grn_models/norman/grnboost2.h5ad \
    --dataset_id norman \
    --ws_consensus resources/grn_benchmark/prior/ws_consensus_norman.csv \
    --ws_distance_background resources/grn_benchmark/prior/ws_distance_background_norman.csv \
    --evaluation_data_sc resources/datasets_raw/norman_sc_counts.h5ad \
    --score output/score.h5ad
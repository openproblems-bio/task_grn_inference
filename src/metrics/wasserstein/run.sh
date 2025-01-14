viash run src/metrics/wasserstein/config.vsh.yaml -- --prediction resources/grn_models/norman/grnboost2.csv \
    --dataset_id norman \
    --ws_consensus resources/prior/ws_consensus_norman.csv \
    --ws_distance_background resources/prior/ws_distance_background_norman.csv \
    --evaluation_data_sc resources/datasets_raw/norman_sc_counts.h5ad \
    --score output/score.h5ad
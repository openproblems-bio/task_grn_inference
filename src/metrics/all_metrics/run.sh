
viash run src/metrics/all_metrics/config.vsh.yaml -- \
    --prediction resources/grn_models/op/scprint.h5ad \
    --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --score output/score.h5ad \
    --regulators_consensus resources/grn_benchmark/prior/regulators_consensus_op.json 
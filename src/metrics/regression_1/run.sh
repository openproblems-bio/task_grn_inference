# viash run src/metrics/regression_1/config.vsh.yaml -- \
#     --prediction resources/grn_models/norman/grnboost2.h5ad \
#     --dataset_id norman --evaluation_data resources/grn_benchmark/evaluation_data//norman_bulk.h5ad \
#     --score output/score.h5ad \
#     --tf_all resources/grn_benchmark/prior/tf_all.csv


python src/metrics/regression_1/script.py \
    --run_local \
    --prediction resources/grn_models/op/grnboost2.h5ad \
    --dataset_id op \
    --evaluation_data resources/grn_benchmark/evaluation_data/op_bulk.h5ad \
    --method_id grnboost2 \
    --score output/score.h5ad 
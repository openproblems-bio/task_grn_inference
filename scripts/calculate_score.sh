# bash src/metrics/all_metrics/run.sh resources/grn_models/norman/grnboost2.csv norman

prediction=${1}
dataset_id=${2}

viash run src/metrics/all_metrics/config.novsh.yaml -- \
    --prediction ${prediction} \
    --dataset_id ${dataset_id} \
    --score output/score.h5ad \
    --tf_all resources/prior/tf_all.csv \
    --regulators_consensus resources/prior/regulators_consensus_${dataset_id}.json \
    --ws_consensus resources/prior/ws_consensus_${dataset_id}.csv \
    --ws_distance_background resources/prior/ws_distance_background_${dataset_id}.csv \
    --evaluation_data_sc resources/evaluation_datasets/${dataset_id}_sc_counts.h5ad \
    --evaluation_data resources/evaluation_datasets/${dataset_id}_perturbation.h5ad
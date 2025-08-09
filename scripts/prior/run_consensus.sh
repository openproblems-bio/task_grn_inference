#!/bin/bash
#SBATCH --job-name=ws
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --time=10:00:00
#SBATCH --qos long
#SBATCH --mem=1000GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

set -e

models_dir="resources/results/test_run/"

echo "Running consensus for regression 2"
# datasets=("op" "norman" "adamson" "nakatake" "replogle" "replogle_scprint" "xaira_HCT116")
datasets=('xaira_HCT116' 'xaira_HEK293T' 'parsebioscience')
models=("pearson_corr" "negative_control" "positive_control" "portia" "ppcor" "scenic" "scprint" "grnboost2" "scenicplus" "scglue" "granie" "figr" "celloracle")

# for dataset in "${datasets[@]}"; do
#     echo "Running consensus for dataset: $dataset"
    
#     python src/metrics/regression_2/consensus/script.py \
#         --dataset "$dataset" \
#         --models_dir "$models_dir" \
#         --regulators_consensus "resources/grn_benchmark/prior/regulators_consensus_${dataset}.json" \
#         --evaluation_data "resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad" \
#         --models "${models[@]}"
# done


# echo "Running consensus for ws distance"

# for dataset in "${datasets[@]}"; do
#     echo "Running consensus for dataset: $dataset"
    
#     python src/metrics/ws_distance/consensus/script.py \
#         --dataset "$dataset" \
#         --models_dir "$models_dir" \
#         --ws_consensus "resources/grn_benchmark/prior/ws_consensus_${dataset}.json" \
#         --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
#         --evaluation_data_sc "resources/grn_benchmark/evaluation_data/${dataset}_sc.h5ad" \
#         --models "${models[@]}"
# done

echo "Calculating scores for all possible connections, WS distance"
for dataset in "${datasets[@]}"; do
    echo "Running consensus for dataset: $dataset"
    
    python src/metrics/ws_distance/background_distance/script.py \
        --dataset "$dataset" \
        --background_distance "resources/grn_benchmark/prior/ws_distance_background_${dataset}.json" \
        --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
        --evaluation_data_sc "resources/grn_benchmark/evaluation_data/${dataset}_sc.h5ad" \
        --max_workers 100 
done

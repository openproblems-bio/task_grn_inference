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
models=("pearson_corr" "negative_control" "positive_control" "portia" "ppcor" "scenic" "scprint" "grnboost2" "scenicplus" "scglue" "granie" "figr" "celloracle")

# datasets=("op" "norman" "adamson" "nakatake" "replogle" "xaira_HCT116" "xaira_HCT116" "xaira_HEK293T" "parsebioscience")

datasets=("parsebioscience")
for dataset in "${datasets[@]}"; do
    echo "Running consensus for dataset: $dataset"
    
    python src/metrics/regression_2/consensus/script.py \
        --dataset "$dataset" \
        --models_dir "$models_dir" \
        --regulators_consensus "resources/grn_benchmark/prior/regulators_consensus_${dataset}.json" \
        --evaluation_data "resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad" \
        --models "${models[@]}"
done


# echo "Running consensus for ws distance"
# # datasets=("op" "norman" "adamson" "nakatake" "replogle" "xaira_HEK293T" "xaira_HCT116" )
# datasets=("xaira_HEK293T" "xaira_HCT116" )
# for dataset in "${datasets[@]}"; do
#     echo "Running consensus for dataset: $dataset"
    
#     python src/metrics/ws_distance/consensus/script.py \
#         --dataset "$dataset" \
#         --models_dir "$models_dir" \
#         --ws_consensus "resources/grn_benchmark/prior/ws_consensus_${dataset}.csv" \
#         --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
#         --evaluation_data_sc "resources/grn_benchmark/evaluation_data/${dataset}_sc.h5ad" \
#         --models "${models[@]}"
# done


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

datasets=("norman" "adamson" "replogle" 'xaira_HEK293T' 'xaira_HCT116' )
# datasets=('xaira_HCT116' )

echo "Calculating scores for all possible connections, WS distance"
for dataset in "${datasets[@]}"; do
    echo "Running for dataset: $dataset"
    
    python src/metrics/ws_distance/background_distance/script.py \
        --dataset "$dataset" \
        --background_distance "resources/grn_benchmark/prior/ws_distance_background_${dataset}.csv" \
        --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
        --evaluation_data_sc "resources/processed_data/${dataset}_evaluation_sc.h5ad" \
        --num_workers 100 
done

#!/bin/bash
#SBATCH --job-name=tf_binding_global
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=30:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

set -euo pipefail

save_dir="output/tf_binding_global"
mkdir -p "$save_dir"

# Global GRN evaluation - only op and 300BCG datasets
datasets=('op' '300BCG')
resources_dir="resources"

# Global GRN files location
global_grn_dir="resources/results/experiment/global_grns"

for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

    # Create separate CSV file for each dataset with ground truth type
    dataset_csv="${save_dir}/summary_${dataset}.csv"
    echo "dataset,method,metric,value" > "$dataset_csv"

    for prediction in ${global_grn_dir}/${dataset}.*.h5ad; do
        if [[ ! -f "$prediction" ]]; then
            echo "  No global GRN files found for dataset $dataset"
            continue
        fi
        
        # Extract method name from filename
        filename=$(basename "$prediction")
        method=$(echo "$filename" | sed "s/${dataset}\.//" | sed 's/\.h5ad$//' | sed 's/\.csv$//')
        
        echo -e "\n  $method  \n"
        mkdir -p "${save_dir}/tmp"
        score="${save_dir}/tmp/${dataset}_${method//[:\/]/_}.h5ad"

        cell_type=$(python -c "
import sys
sys.path.insert(0, 'src/utils')
from config import DATASETS_CELLTYPES
print(DATASETS_CELLTYPES.get('$dataset', ''))
")
        python src/metrics/tf_binding/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
            --ground_truth_unibind ${resources_dir}/grn_benchmark/ground_truth/${cell_type}_unibind.csv \
            --ground_truth_chipatlas ${resources_dir}/grn_benchmark/ground_truth/${cell_type}_chipatlas.csv \
            --ground_truth_remap ${resources_dir}/grn_benchmark/ground_truth/${cell_type}_remap.csv \
            --score "$score"

    done  # end global GRN files loop

    echo -e "\nResults for dataset $dataset collected in: $dataset_csv"
done  # end datasets loop

echo -e "\nAll dataset results saved in separate CSV files in: $save_dir"

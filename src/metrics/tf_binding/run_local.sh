#!/bin/bash
#SBATCH --job-name=tf_binding
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

save_dir="output/tf_binding"
mkdir -p "$save_dir"

# datasets to process
datasets=(  'replogle' 'norman'  'adamson'  "300BCG" "ibd" 'parsebioscience' 'op' ) #'xaira_HCT116'
# datasets=( 'xaira_HEK293T') 
resources_dir="resources"
# methods to process - focus on granie for debugging

methods=(  "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie" "scgpt" "geneformer")
# methods=( "pearson_corr" "negative_control" "positive_control" )
for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

    # Create separate CSV file for each dataset with ground truth type
    dataset_csv="${save_dir}/summary_${dataset}.csv"
    echo "dataset,method,metric,value" > "$dataset_csv"

    for method in "${methods[@]}"; do
        echo -e "\n  $method  \n"
        prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"
        mkdir -p "${save_dir}/tmp"
        score="${save_dir}/tmp/${dataset}_${method}.h5ad"

        if [[ ! -f "$prediction" ]]; then
            echo "File not found: $prediction, skipping..."
            continue
        fi

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


    done  # end methods loop

    echo -e "\nResults for dataset $dataset collected in: $dataset_csv"
done  # end datasets loop

echo -e "\nAll dataset results saved in separate CSV files in: $save_dir"
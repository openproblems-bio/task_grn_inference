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
datasets=(  'replogle' 'norman'  'adamson'  "300BCG" 'ibd_uc' 'ibd_cd' 'parsebioscience' 'op' 'xaira_HCT116' 'xaira_HEK293T' ) #'xaira_HCT116'
# datasets=(  'replogle'  ) #'xaira_HCT116'
resources_dir="resources"
methods=(  "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie" "scgpt" "geneformer")
# methods=( "grnboost" )


# Create summary CSV file with header - will be written by first Python call
summary_csv="${save_dir}/summary.csv"
# Header will be dynamically created from metric names
echo -n "" > "$summary_csv"
first_run=true

for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"


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

        # Extract metrics from the .h5ad and append to CSV
        python -u - <<EOF
import anndata as ad
import pandas as pd
import numpy as np
import os

adata = ad.read_h5ad("${score}")
if "metric_values" in adata.uns:
    metric_names = adata.uns["metric_ids"]
    metric_values = adata.uns["metric_values"]
    
    # Single row format - transpose the metrics into columns
    df = pd.DataFrame([metric_values], columns=metric_names)
    df.insert(0, "dataset", "${dataset}")
    df.insert(1, "method", "${method}")
    
    # Write header if file is empty
    write_header = not os.path.exists("${summary_csv}") or os.path.getsize("${summary_csv}") == 0
    df.to_csv("${summary_csv}", mode="a", header=write_header, index=False)
EOF

    done  # end methods loop

done  # end datasets loop

echo -e "\nAll results saved in: $summary_csv"
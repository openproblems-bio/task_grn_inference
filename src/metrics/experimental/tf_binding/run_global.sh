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

save_dir="resources/results/experiment/tf_binding"
mkdir -p "$save_dir"

# Global GRN evaluation - only op and 300BCG datasets
datasets=('op')
resources_dir="resources"

# Global GRN files location
global_grn_dir="resources/results/experiment/global_grns"

# Create summary CSV file
summary_csv="${save_dir}/summary_global.csv"
echo "dataset,method,gt_source,metric,value" > "$summary_csv"

for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

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

        python src/metrics/experimental/tf_binding/script.py \
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

adata = ad.read_h5ad("${score}")
if "metric_values" in adata.uns:
    metric_names = adata.uns["metric_ids"]
    metric_values = adata.uns["metric_values"]
    gt_sources = adata.uns.get("gt_sources", [])
    
    # Check if metric_values is 2D array (multiple rows) or 1D (single row)
    if isinstance(metric_values, np.ndarray) and metric_values.ndim == 2:
        # Multiple rows (one per ground truth source)
        for row_idx, row_values in enumerate(metric_values):
            df = pd.DataFrame({"metric": metric_names, "value": row_values})
            df["dataset"] = "${dataset}"
            method_clean = "${method}".replace(',', ';')
            df["method"] = method_clean
            df["gt_source"] = gt_sources[row_idx] if row_idx < len(gt_sources) else ""
            df = df[["dataset", "method", "gt_source", "metric", "value"]]
            df.to_csv("${summary_csv}", mode="a", header=False, index=False)
    else:
        # Single row (legacy format)
        df = pd.DataFrame({"metric": metric_names, "value": metric_values})
        df["dataset"] = "${dataset}"
        method_clean = "${method}".replace(',', ';')
        df["method"] = method_clean
        df["gt_source"] = gt_sources[0] if len(gt_sources) > 0 else ""
        df = df[["dataset", "method", "gt_source", "metric", "value"]]
        df.to_csv("${summary_csv}", mode="a", header=False, index=False)
EOF

    done  # end global GRN files loop

done  # end datasets loop

echo -e "\nAll results saved in: $summary_csv"

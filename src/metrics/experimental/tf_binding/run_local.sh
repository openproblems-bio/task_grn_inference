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
datasets=(  'replogle' 'norman'  'adamson'  "300BCG" "ibd" 'parsebioscience' 'op' 'xaira_HCT116' 'xaira_HEK293T' ) #'xaira_HCT116'
# datasets=(  'replogle'  ) #'xaira_HCT116'
resources_dir="resources"
methods=(  "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie" "scgpt" "geneformer")
# methods=( "grnboost" )


# Create summary CSV file
summary_csv="${save_dir}/summary.csv"
echo "dataset,method,gt_source,metric,value" > "$summary_csv"

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
            df["method"] = "${method}"
            df["gt_source"] = gt_sources[row_idx] if row_idx < len(gt_sources) else ""
            df = df[["dataset", "method", "gt_source", "metric", "value"]]
            df.to_csv("${summary_csv}", mode="a", header=False, index=False)
    else:
        # Single row (legacy format)
        df = pd.DataFrame({"metric": metric_names, "value": metric_values})
        df["dataset"] = "${dataset}"
        df["method"] = "${method}"
        df["gt_source"] = gt_sources[0] if len(gt_sources) > 0 else ""
        df = df[["dataset", "method", "gt_source", "metric", "value"]]
        df.to_csv("${summary_csv}", mode="a", header=False, index=False)
EOF

    done  # end methods loop

done  # end datasets loop

echo -e "\nAll results saved in: $summary_csv"
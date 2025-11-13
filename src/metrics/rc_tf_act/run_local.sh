#!/bin/bash
#SBATCH --job-name=rc_tf_act
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

save_dir="output/rc_tf_act"
mkdir -p "$save_dir"

# Datasets to process (only those with rc_tf_ac grouping defined)
datasets=('300BCG'  'parsebioscience'  'op')
# datasets=('op' )
# Methods to process
methods=("grnboost" "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie")
# methods=("grnboost")
# Combined CSV for all results
combined_csv="${save_dir}/rc_tf_act_scores.csv"
echo "dataset,method,metric,value" > "$combined_csv"

for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"
    
    # Create separate CSV file for each dataset
    dataset_csv="${save_dir}/rc_tf_act_scores_${dataset}.csv"
    echo "dataset,method,metric,value" > "$dataset_csv"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

    for method in "${methods[@]}"; do
        prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"
        score="${save_dir}/rc_tf_act_${dataset}_${method}.h5ad"

        if [[ ! -f "$prediction" ]]; then
            echo "File not found: $prediction, skipping..."
            continue
        fi

        echo -e "\nProcessing method: $method\n"
        python src/metrics/rc_tf_act/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
            --score "$score" \
            --min_targets 5
        
        # Extract metrics from the .h5ad and append to CSV
        python -u - <<EOF
import anndata as ad
import pandas as pd

adata = ad.read_h5ad("${score}")
if "metric_values" in adata.uns:
    metric_names = adata.uns["metric_ids"]
    metric_values = adata.uns["metric_values"]
    df = pd.DataFrame({"metric": metric_names, "value": metric_values})
    df["dataset"] = "${dataset}"
    df["method"] = "${method}"
    df = df[["dataset", "method", "metric", "value"]]  # Reorder columns
    df.to_csv("${dataset_csv}", mode="a", header=False, index=False)
    df.to_csv("${combined_csv}", mode="a", header=False, index=False)
EOF

    done
    echo -e "\nResults for dataset $dataset collected in: $dataset_csv"
done

echo -e "\nAll results saved in: $combined_csv"

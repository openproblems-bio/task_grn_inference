#!/bin/bash
#SBATCH --job-name=ar_global
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

datasets=('op' '300BCG')
save_dir='output/ar_global/'
mkdir -p "$save_dir"

# Global GRN files location
global_grn_dir="resources/results/experiment/global_grns"

for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"
    dataset_csv="${save_dir}/ar_scores_${dataset}.csv"
    echo "dataset,method,metric,value" > "$dataset_csv"
    
    for prediction in ${global_grn_dir}/${dataset}.*.h5ad; do
        if [[ ! -f "$prediction" ]]; then
            echo "  No global GRN files found for dataset $dataset"
            continue
        fi
        
        filename=$(basename "$prediction")
        method=$(echo "$filename" | sed "s/${dataset}\.//" | sed 's/\.h5ad$//' | sed 's/\.csv$//')
        
        score="${save_dir}/ar_${dataset}_${method//[:\/]/_}.h5ad"

        echo -e "\nProcessing method: $method\n"
        python src/metrics/ar/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
            --score "$score"

        python -u - <<EOF
import anndata as ad
import pandas as pd

adata = ad.read_h5ad("${score}")
if "metric_values" in adata.uns:
    metric_names = adata.uns["metric_ids"]
    metric_values = adata.uns["metric_values"]
    df = pd.DataFrame({"metric": metric_names, "value": metric_values})
    df["dataset"] = "${dataset}"
    method_clean = "${method}".replace(',', ';')
    df["method"] = method_clean
    df = df[["dataset", "method", "metric", "value"]]
    df.to_csv("${dataset_csv}", mode="a", header=False, index=False)
EOF

    done
    echo -e "\nResults for dataset $dataset collected in: $dataset_csv"
done

echo -e "\nAll dataset results saved in separate CSV files in: $save_dir"

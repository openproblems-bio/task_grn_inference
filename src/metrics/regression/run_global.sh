#!/bin/bash
#SBATCH --job-name=reg2_global
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

set -euo pipefail

save_dir="output/regression"
mkdir -p "$save_dir"

# Global GRN evaluation - only op and 300BCG datasets
datasets=('op' '300BCG')

# Global GRN files location
global_grn_dir="resources/results/experiment/global_grns"

# Create summary CSV file
summary_csv="${save_dir}/summary_global.csv"
echo "dataset,method,metric,value" > "$summary_csv"

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
        
        score="${save_dir}/reg2_${dataset}_${method//[:\/]/_}.h5ad"

        echo -e "\nProcessing method: $method\n"
        python src/metrics/regression/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
            --regulators_consensus "resources/grn_benchmark/prior/regulators_consensus_${dataset}.json" \
            --score "$score"
        
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
    method_clean = "${method}".replace(',', ';')
    df["method"] = method_clean
    df = df[["dataset", "method", "metric", "value"]]
    df.to_csv("${summary_csv}", mode="a", header=False, index=False)
EOF

    done
done

echo -e "\nAll results saved in: $summary_csv"

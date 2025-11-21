#!/bin/bash
#SBATCH --job-name=rc
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

# datasets=( "ibd" "parsebioscience" "300BCG" "op")
datasets=('op') 
methods=( "pearson_corr" "negative_control" "positive_control"  "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie")
# methods=("grnboost")

# temporary file to collect CSV rows
save_dir='output/rc/'
mkdir -p "$save_dir"

# Create summary CSV file
summary_csv="${save_dir}/summary.csv"
echo "dataset,method,metric,value" > "$summary_csv"

for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"
    for method in "${methods[@]}"; do
        prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"
        score="${save_dir}/sem_${dataset}_${method}.h5ad"

        if [[ ! -f "$prediction" ]]; then
            echo "File not found: $prediction, skipping..."
            continue
        fi

        echo -e "\nProcessing method: $method\n"
        python src/metrics/rc/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
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
    df["method"] = "${method}"
    df = df[["dataset", "method", "metric", "value"]]  # Reorder columns to match header
    df.to_csv("${summary_csv}", mode="a", header=False, index=False)
EOF

    done
done

echo -e "\nAll results saved in: $summary_csv"
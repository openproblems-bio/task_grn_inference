#!/bin/bash
#SBATCH --job-name=recovery_2
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

datasets=( 'op') #"300BCG" "ibd" 'parsebioscience'
methods=("negative_control" "pearson_corr" "positive_control"  "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie")
# methods=("granie")

# temporary file to collect CSV rows
save_dir='output/temp/'
combined_csv="${save_dir}/scores.csv"
echo "dataset,method,metric,value" > "$combined_csv"

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
        python src/metrics/recovery_2/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
            --score "$score"

        # Extract metrics from the .h5ad and append to CSV
        python - <<EOF
import anndata as ad
import pandas as pd

adata = ad.read_h5ad("${score}")
if "metrics_value" in adata.uns:
    metrics = adata.uns["metrics_value"]
    df = pd.DataFrame(list(metrics.items()), columns=["metric", "value"])
    df["dataset"] = "${dataset}"
    df["method"] = "${method}"
    df = df[["dataset", "method", "metric", "value"]]
    df.to_csv("${combined_csv}", mode="a", header=False, index=False)
EOF

    done
done

echo -e "\nAll results collected in: $combined_csv"
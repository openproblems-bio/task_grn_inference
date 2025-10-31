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
datasets=(  'replogle' 'norman'  'adamson' ) # 'xaira_HCT116' 'replogle' 'norman'  'adamson') #"300BCG" "ibd" 'parsebioscience''op' "300BCG" 'parsebioscience'   'replogle' 'norman' 'adamson'
# methods to process
methods=( "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie")

for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"
    
    # Create separate CSV file for each dataset
    dataset_csv="${save_dir}/tf_binding_scores_${dataset}.csv"
    echo "dataset,method,metric,value" > "$dataset_csv"

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

    for method in "${methods[@]}"; do
        prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"
        score="${save_dir}/tf_binding_${dataset}_${method}.h5ad"

        if [[ ! -f "$prediction" ]]; then
            echo "File not found: $prediction, skipping..."
            continue
        fi
        if [[ "$dataset" == "replogle" || "$dataset" == "norman" || "$dataset" == "adamson" ]]; then
            ground_truth="resources/grn_benchmark/ground_truth/K562_remap.csv"
        elif [[ "$dataset" == "xaira_HEK293T" ]]; then
            ground_truth="resources/grn_benchmark/ground_truth/HEK293T_remap.csv"
        elif [[ "$dataset" == "xaira_HCT116" ]]; then
            ground_truth="resources/grn_benchmark/ground_truth/HCT116_chipatlas.csv"
        else
            echo "No ground truth available for dataset: $dataset, skipping..."
            continue
        fi

        echo -e "\nProcessing method: $method\n"
        python src/metrics/tf_binding/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
            --ground_truth "$ground_truth" \
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
    df.to_csv("${dataset_csv}", mode="a", header=False, index=False)
EOF

    done
    echo -e "\nResults for dataset $dataset collected in: $dataset_csv"
done

echo -e "\nAll dataset results saved in separate CSV files in: $save_dir"
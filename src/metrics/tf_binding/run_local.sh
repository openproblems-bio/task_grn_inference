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
# datasets=(  'replogle' 'norman'  'adamson'  "300BCG" "ibd" 'parsebioscience' 'op' ) #'xaira_HCT116'
datasets=( 'xaira_HEK293T' 'xaira_HCT116' 'replogle' 'norman'  'adamson' 'op'  "300BCG" "ibd" 'parsebioscience' ) 

# methods to process - focus on granie for debugging

methods=(  "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie")
# methods=( "pearson_corr" "negative_control" "positive_control" )
for dataset in "${datasets[@]}"; do
    echo -e "\n\nProcessing dataset: $dataset\n"
    
    

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

    # Define ground truth types to process
    ground_truth_types=("remap2022" "chipatlas" "unibind")

    for gt_type in "${ground_truth_types[@]}"; do
        echo -e "\n--- $gt_type ---\n"

        # Create separate CSV file for each dataset with ground truth type
        dataset_csv="${save_dir}/summary_${dataset}_${gt_type}.csv"
        echo "dataset,method,ground_truth_type,metric,value" > "$dataset_csv"

        for method in "${methods[@]}"; do
            echo -e "\n  $method  \n"
            prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"
            mkdir -p "${save_dir}/tmp"
            score="${save_dir}/tmp/${dataset}_${method}_${gt_type}.h5ad"

            if [[ ! -f "$prediction" ]]; then
                echo "File not found: $prediction, skipping..."
                continue
            fi

            # Set ground truth file based on dataset and GT type
            if [[ "$dataset" == "replogle" || "$dataset" == "norman" || "$dataset" == "adamson" ]]; then
                # K562 cell line datasets
                ground_truth="resources/grn_benchmark/ground_truth/K562_${gt_type}_k562.csv"
            elif [[ "$dataset" == "xaira_HEK293T" ]]; then
                # HEK293 datasets  
                ground_truth="resources/grn_benchmark/ground_truth/HEK293_${gt_type}_hek293.csv"
            elif [[ "$dataset" == "xaira_HCT116" ]]; then
                # HCT116 datasets
                ground_truth="resources/grn_benchmark/ground_truth/HCT116_${gt_type}_hct116.csv"
            elif [[ "$dataset" == "op" || "$dataset" == "parsebioscience" || "$dataset" == "300BCG" || "$dataset" == "ibd" ]]; then
                # PBMC datasets
                ground_truth="resources/grn_benchmark/ground_truth/PBMC_${gt_type}_pbmc.csv"
            else
                echo "No ground truth mapping defined for dataset: $dataset, skipping..."
                continue
            fi

            # Check if ground truth file exists
            if [[ ! -f "$ground_truth" ]]; then
                echo "Ground truth file not found: $ground_truth, skipping..."
                continue
            fi

            echo "Using ground truth: $ground_truth"
            
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
    df["ground_truth_type"] = "${gt_type}"
    df = df[["dataset", "method", "ground_truth_type", "metric", "value"]]  # Reorder columns to match header
    df.to_csv("${dataset_csv}", mode="a", header=False, index=False)
    
    # Print summary for this method and GT type
    precision = df[df['metric'] == 'tf_binding_precision']['value'].iloc[0] if 'tf_binding_precision' in df['metric'].values else 0
    recall = df[df['metric'] == 'tf_binding_recall']['value'].iloc[0] if 'tf_binding_recall' in df['metric'].values else 0
    n_tfs_evaluated = df[df['metric'] == 'n_tfs_evaluated']['value'].iloc[0] if 'n_tfs_evaluated' in df['metric'].values else 0
    n_tfs_in_gt = df[df['metric'] == 'n_tfs_in_gt']['value'].iloc[0] if 'n_tfs_in_gt' in df['metric'].values else 0
    tf_coverage = df[df['metric'] == 'tf_coverage_pct']['value'].iloc[0] if 'tf_coverage_pct' in df['metric'].values else 0
    
    # print(f"\n=== ${method} - ${gt_type} SUMMARY ===")
    # print(f"Precision: {precision:.4f}")
    # print(f"Recall: {recall:.4f}")
    # print(f"TFs evaluated: {n_tfs_evaluated}/{n_tfs_in_gt} ({tf_coverage:.1f}%)")
else:
    print("No metrics found in ${score}")
EOF

        done  # end methods loop
    done  # end ground truth types loop
    echo -e "\nResults for dataset $dataset collected in: $dataset_csv"
done  # end datasets loop

echo -e "\nAll dataset results saved in separate CSV files in: $save_dir"
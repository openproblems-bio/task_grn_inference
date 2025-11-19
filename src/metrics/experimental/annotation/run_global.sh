#!/bin/bash
#SBATCH --job-name=pathway_annotation_global
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

save_dir="output/pathway_annotation_global"
mkdir -p "$save_dir"

# Pathway file location
pathway_file="/vol/projects/jnourisa/prior/h.all.v2024.1.Hs.symbols.gmt"

# datasets to process (only op and 300BCG for global GRN evaluation)
datasets=( 'op' '300BCG' )

resources_dir="resources"

# Global GRN files location
global_grn_dir="resources/results/experiment/global_grns"

for dataset in "${datasets[@]}"; do
    echo -e "\n\n=========================================="
    echo "Processing dataset: $dataset"
    echo "==========================================" 

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

    # Create separate CSV file for each dataset
    dataset_csv="${save_dir}/pathway_annotation_scores_${dataset}.csv"
    echo "dataset,method,metric,value" > "$dataset_csv"

    # Process all global GRN files for this dataset
    for prediction_file in ${global_grn_dir}/${dataset}.*.h5ad; do
        if [[ ! -f "$prediction_file" ]]; then
            echo "  No global GRN files found for dataset $dataset"
            continue
        fi
        
        # Extract method name from filename
        # Format: op.METHOD.h5ad or op.METHOD:TISSUE.csv.h5ad
        filename=$(basename "$prediction_file")
        method=$(echo "$filename" | sed "s/${dataset}\.//" | sed 's/\.h5ad$//' | sed 's/\.csv$//')
        
        echo -e "\n  Processing method: $method"
        
        mkdir -p "${save_dir}/tmp"
        score="${save_dir}/tmp/${dataset}_${method//[:\/]/_}.h5ad"

        # Run pathway annotation metric
        python src/metrics/experimental/annotation/script.py \
            --prediction "$prediction_file" \
            --evaluation_data "$evaluation_data" \
            --tf_all ${resources_dir}/grn_benchmark/prior/tf_all.csv \
            --pathway_file "$pathway_file" \
            --score "$score" \
            --fdr_threshold 0.05 \
            --min_targets 10 \
            --max_targets 100 \
            --run_gene_set_recovery \
            --ulm_activity_threshold 0.0 \
            --ulm_pvalue_threshold 0.01 \
            --ulm_baseline_method zero_centered

        # Read the scores from the h5ad file and append to CSV
        if [[ -f "$score" ]]; then
            python -c "
import anndata as ad
import sys

try:
    adata = ad.read_h5ad('$score')
    metric_ids = adata.uns.get('metric_ids', [])
    metric_values = adata.uns.get('metric_values', [])
    
    method_clean = '$method'.replace(',', ';')  # Replace commas in method name
    
    for metric_id, metric_value in zip(metric_ids, metric_values):
        print(f'$dataset,{method_clean},{metric_id},{metric_value}')
except Exception as e:
    print(f'Error reading scores: {e}', file=sys.stderr)
    sys.exit(1)
" >> "$dataset_csv"
            echo "    Scores appended to CSV"
        else
            echo "    Warning: Score file not created"
        fi

    done  # end global GRN files loop

    echo -e "\nResults for dataset $dataset saved to: $dataset_csv"
done  # end datasets loop

echo -e "\n=========================================="
echo "All results saved in: $save_dir"
echo "=========================================="

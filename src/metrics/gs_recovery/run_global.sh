#!/bin/bash
#SBATCH --job-name=gs_recovery_global
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

save_dir="resources/results/experiment/gs_recovery"
mkdir -p "$save_dir"

# Pathway files location
pathway_dir="resources/grn_benchmark/prior/pathways"

datasets=( 'op' )

resources_dir="resources"

# Global GRN files location
global_grn_dir="resources/results/experiment/global_grns"

# Create summary CSV file
summary_csv="${save_dir}/summary_global.csv"
echo "dataset,method,metric,value" > "$summary_csv"

for dataset in "${datasets[@]}"; do
    echo -e "\n\n=========================================="
    echo "Processing dataset: $dataset"
    echo "==========================================" 

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

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

        # Build pathway files argument (comma-separated list of geneset_name:file_path)
        pathway_files_arg=""
        for geneset_file in ${pathway_dir}/*.gmt; do
            if [[ -f "$geneset_file" ]]; then
                geneset_name=$(basename "$geneset_file" .gmt)
                if [[ -z "$pathway_files_arg" ]]; then
                    pathway_files_arg="${geneset_name}:${geneset_file}"
                else
                    pathway_files_arg="${pathway_files_arg},${geneset_name}:${geneset_file}"
                fi
            fi
        done

        # Run pathway gs_recovery metric
        python src/metrics/experimental/gs_recovery/script.py \
            --prediction "$prediction_file" \
            --evaluation_data "$evaluation_data" \
            --tf_all ${resources_dir}/grn_benchmark/prior/tf_all.csv \
            --pathway_files "$pathway_files_arg" \
            --score "$score" \
            --fdr_threshold 0.05 \
            --min_targets 10 \
            --max_targets 100 \
            --run_gene_set_recovery \
            --ulm_activity_threshold 0.0 \
            --ulm_pvalue_threshold 0.01 \
            --ulm_baseline_method zero_centered \
            --n_workers 20

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
" >> "$summary_csv"
            echo "    Scores appended to CSV"
        else
            echo "    Warning: Score file not created"
        fi

    done  # end global GRN files loop

done  # end datasets loop

echo -e "\n=========================================="
echo "All results saved in: $summary_csv"
echo "=========================================="

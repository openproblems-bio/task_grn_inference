#!/bin/bash
#SBATCH --job-name=pathway_annotation
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

save_dir="output/annotation"
mkdir -p "$save_dir"

# Pathway file location
pathway_file="/vol/projects/jnourisa/prior/h.all.v2024.1.Hs.symbols.gmt"

# datasets to process
datasets=( 'op' 'replogle' 'norman' 'adamson' "300BCG" "ibd" 'parsebioscience' 'xaira_HCT116' 'xaira_HEK293T' )
# datasets=( 'op' )  # Test with single dataset first

resources_dir="resources"
methods=( "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie" "scgpt" "geneformer" )
# methods=( "pearson_corr" "grnboost" )  # Test subset

# Create summary CSV file
summary_csv="${save_dir}/summary.csv"
echo "dataset,method,metric,value" > "$summary_csv"

for dataset in "${datasets[@]}"; do
    echo -e "\n\n=========================================="
    echo "Processing dataset: $dataset"
    echo "==========================================" 

    evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

    for method in "${methods[@]}"; do
        echo -e "\n  Processing method: $method"
        prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"
        mkdir -p "${save_dir}/tmp"
        score="${save_dir}/tmp/${dataset}_${method}.h5ad"

        if [[ ! -f "$prediction" ]]; then
            echo "    File not found: $prediction, skipping..."
            continue
        fi

        # Run pathway annotation metric
        python src/metrics/experimental/annotation/script.py \
            --prediction "$prediction" \
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
    
    for metric_id, metric_value in zip(metric_ids, metric_values):
        print(f'$dataset,$method,{metric_id},{metric_value}')
except Exception as e:
    print(f'Error reading scores: {e}', file=sys.stderr)
    sys.exit(1)
" >> "$summary_csv"
            echo "    Scores appended to CSV"
        else
            echo "    Warning: Score file not created"
        fi

    done  # end methods loop

done  # end datasets loop

echo -e "\n=========================================="
echo "All results saved in: $summary_csv"
echo "=========================================="
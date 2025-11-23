#!/bin/bash
#SBATCH --job-name=pathway_gs_recovery
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

save_dir="output/gs_recovery"
mkdir -p "$save_dir"

# Pathway files - define individual paths
resources_dir="resources"
pathway_dir="${resources_dir}/grn_benchmark/prior/pathways"

geneset_hallmark_2020="${pathway_dir}/hallmark_2020.gmt"
geneset_kegg_2021="${pathway_dir}/kegg_2021.gmt"
geneset_reactome_2022="${pathway_dir}/reactome_2022.gmt"
geneset_go_bp_2023="${pathway_dir}/go_bp_2023.gmt"
geneset_bioplanet_2019="${pathway_dir}/bioplanet_2019.gmt"
geneset_wikipathways_2019="${pathway_dir}/wikipathways_2019.gmt"

# datasets to process
datasets=( 'op' 'replogle' 'norman' 'adamson' "300BCG" "ibd" 'parsebioscience' 'xaira_HCT116' 'xaira_HEK293T' )
# datasets to process
datasets=( 'op' 'replogle' 'norman' 'adamson' "300BCG" "ibd" 'parsebioscience' 'xaira_HCT116' 'xaira_HEK293T' )
# datasets=( 'op' )  # Test with single dataset first

methods=( "pearson_corr" "negative_control" "positive_control" "ppcor" "portia" "scenic" "grnboost" "scprint" "scenicplus" "celloracle" "scglue" "figr" "granie" "scgpt" "geneformer" )
# methods=( "pearson_corr" "grnboost" )  # Test subset

# Create summary CSV file with dynamic header
summary_csv="${save_dir}/summary.csv"
echo -n "" > "$summary_csv"  # Will write header dynamically from first result

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

        # Run pathway gs_recovery metric
        python src/metrics/gs_recovery/script.py \
            --prediction "$prediction" \
            --evaluation_data "$evaluation_data" \
            --tf_all ${resources_dir}/grn_benchmark/prior/tf_all.csv \
            --geneset_hallmark_2020 "$geneset_hallmark_2020" \
            --geneset_kegg_2021 "$geneset_kegg_2021" \
            --geneset_reactome_2022 "$geneset_reactome_2022" \
            --geneset_go_bp_2023 "$geneset_go_bp_2023" \
            --geneset_bioplanet_2019 "$geneset_bioplanet_2019" \
            --geneset_wikipathways_2019 "$geneset_wikipathways_2019" \
            --score "$score" 

        # Read the scores from the h5ad file and append to CSV (single-row format)
        if [[ -f "$score" ]]; then
            python -c "
import anndata as ad
import pandas as pd
import sys
import os

try:
    adata = ad.read_h5ad('$score')
    metric_ids = adata.uns.get('metric_ids', [])
    metric_values = adata.uns.get('metric_values', [])
    
    # Create single-row DataFrame
    df = pd.DataFrame([metric_values], columns=metric_ids)
    df.insert(0, 'dataset', '$dataset')
    df.insert(1, 'method', '$method')
    
    # Write header if file is empty
    write_header = not os.path.exists('$summary_csv') or os.path.getsize('$summary_csv') == 0
    df.to_csv('$summary_csv', mode='a', header=write_header, index=False)
    
except Exception as e:
    print(f'Error reading scores: {e}', file=sys.stderr)
    sys.exit(1)
" 
            echo "    Scores appended to CSV"
        else
            echo "    Warning: Score file not created"
        fi

    done  # end methods loop

done  # end datasets loop

echo -e "\n=========================================="
echo "All results saved in: $summary_csv"
echo "=========================================="
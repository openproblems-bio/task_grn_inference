#!/bin/bash
# Test script to compare original vs improved VC implementation on real data

set -euo pipefail

# Create output directories
mkdir -p "output/vc_comparison/original"
mkdir -p "output/vc_comparison/improved"

# Test on a subset of methods for quick comparison
dataset="op"
test_methods=("positive_control" "negative_control" "pearson_corr" "grnboost")

evaluation_data="resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad"

echo "=== Testing VC Implementations on Real Data ==="
echo "Dataset: $dataset"
echo "Methods: ${test_methods[*]}"
echo ""

# CSV files to collect results
original_csv="output/vc_comparison/original_vc_scores.csv"
improved_csv="output/vc_comparison/improved_vc_scores.csv"

echo "method,vc_score" > "$original_csv"
echo "method,vc_score" > "$improved_csv"

for method in "${test_methods[@]}"; do
    prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"
    
    if [[ ! -f "$prediction" ]]; then
        echo "File not found: $prediction, skipping..."
        continue
    fi

    echo "Processing method: $method"

    # Run original VC implementation
    echo "  Running original VC..."
    original_score="output/vc_comparison/original/vc_${dataset}_${method}.h5ad"
    python src/metrics/vc/script.py \
        --prediction "$prediction" \
        --evaluation_data "$evaluation_data" \
        --score "$original_score" 2>/dev/null || echo "Original VC failed for $method"

    # Run improved VC implementation  
    echo "  Running improved VC..."
    improved_score="output/vc_comparison/improved/vc_${dataset}_${method}.h5ad"
    python src/metrics/vc/improved_script.py \
        --prediction "$prediction" \
        --evaluation_data "$evaluation_data" \
        --score "$improved_score" \
        --use_improved true 2>/dev/null || echo "Improved VC failed for $method"

    echo "  Completed $method"
    echo ""
done

echo "=== Comparison Complete ==="
echo "Results saved in output/vc_comparison/"
echo ""

# Extract and compare scores
echo "=== Score Comparison ==="
printf "%-20s %-15s %-15s %-15s\n" "Method" "Original VC" "Improved VC" "Difference"
printf "%-20s %-15s %-15s %-15s\n" "--------------------" "---------------" "---------------" "---------------"

for method in "${test_methods[@]}"; do
    original_score_file="output/vc_comparison/original/vc_${dataset}_${method}.h5ad"
    improved_score_file="output/vc_comparison/improved/vc_${dataset}_${method}.h5ad"
    
    if [[ -f "$original_score_file" ]] && [[ -f "$improved_score_file" ]]; then
        # Extract scores from h5ad files using Python
        original_score=$(python -c "
import anndata as ad
import sys
try:
    adata = ad.read_h5ad('$original_score_file')
    score = float(adata.uns['metric_meta']['vc'])
    print(f'{score:.6f}')
except:
    print('ERROR')
")
        
        improved_score=$(python -c "
import anndata as ad
import sys
try:
    adata = ad.read_h5ad('$improved_score_file')
    score = float(adata.uns['metric_meta']['vc'])
    print(f'{score:.6f}')
except:
    print('ERROR')
")
        
        if [[ "$original_score" != "ERROR" ]] && [[ "$improved_score" != "ERROR" ]]; then
            difference=$(python -c "print(f'{float('$improved_score') - float('$original_score'):.6f}')")
            printf "%-20s %-15s %-15s %-15s\n" "$method" "$original_score" "$improved_score" "$difference"
            
            # Add to CSV files
            echo "$method,$original_score" >> "$original_csv"
            echo "$method,$improved_score" >> "$improved_csv"
        else
            printf "%-20s %-15s %-15s %-15s\n" "$method" "ERROR" "ERROR" "N/A"
        fi
    else
        printf "%-20s %-15s %-15s %-15s\n" "$method" "MISSING" "MISSING" "N/A"
    fi
done

echo ""
echo "Detailed results saved in:"
echo "  Original: $original_csv"
echo "  Improved: $improved_csv"
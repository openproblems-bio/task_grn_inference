#!/bin/bash

# Function to parse command line arguments
parse_arguments() {
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --rna)
                rna="$2"
                shift 2
                ;;
            --prediction)
                prediction="$2"
                shift 2
                ;;
            --atac)
                atac="$2"
                shift 2
                ;;
            --rna_all)
                rna_all="$2"
                shift 2
                ;;
            --layer)
                layer="$2"
                shift 2
                ;;
            *)  # Ignore unknown arguments
                shift
                ;;
        esac
    done

    # Check if required arguments are provided (only when not using dataset)
    if [[ -z "$prediction" ]]; then
        echo "Error: --prediction argument is required"
        exit 1
    fi

    # For non-ATAC methods, at least one of rna or rna_all should be provided
    if [[ -z "$rna" ]] && [[ -z "$rna_all" ]]; then
        echo "Warning: Neither --rna nor --rna_all provided. Some methods may require one of these."
    fi
    
}
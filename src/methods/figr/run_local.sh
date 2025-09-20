#!/bin/bash
#SBATCH --job-name=figr
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=500GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

method="figr"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments
parse_arguments "$@"

# Pass arguments to R script
r_args="--rna $rna --prediction $prediction"
if [ ! -z "$atac" ]; then
    r_args="$r_args --atac $atac"
fi
if [ ! -z "$layer" ]; then
    r_args="$r_args --layer $layer"
fi

singularity run ../../images/${method} Rscript src/methods/${method}/script.R $r_args
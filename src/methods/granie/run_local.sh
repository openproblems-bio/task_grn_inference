#!/bin/bash
#SBATCH --job-name=granie
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

method="granie"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments
parse_arguments "$@"

# Derive dataset name from prediction path for dataset-specific temp dir
dataset=$(basename $(dirname "$prediction"))
granie_temp_dir="output/granie/${dataset}"

# Pass arguments to R script
r_args="--rna $rna --atac $atac --prediction $prediction --temp_dir $granie_temp_dir"
if [ ! -z "$layer" ]; then
    r_args="$r_args --layer $layer"
fi

export SINGULARITYENV_RETICULATE_PYTHON=/usr/bin/python3
singularity run resources/singularity/${method} Rscript src/methods/${method}/script.R $r_args
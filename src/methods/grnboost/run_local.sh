#!/bin/bash
#SBATCH --job-name=grnboost2
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=25:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

method="grnboost"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments
parse_arguments "$@"

# Pass arguments to Python script
python_args="--rna $rna --prediction $prediction"
if [ ! -z "$layer" ]; then
    python_args="$python_args --layer $layer"
fi

singularity run ../../images/scenic python src/methods/${method}/script.py $python_args
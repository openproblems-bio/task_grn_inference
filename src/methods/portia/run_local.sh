#!/bin/bash
#SBATCH --job-name=portia
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

method="portia"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments
parse_arguments "$@"

# Pass arguments to Python script
python_args="--prediction $prediction"
if [ ! -z "$rna" ]; then
    python_args="$python_args --rna $rna"
fi
if [ ! -z "$layer" ]; then
    python_args="$python_args --layer $layer"
fi

singularity run ../../images/${method} python src/methods/${method}/script.py $python_args
#!/bin/bash
#SBATCH --job-name=scglue
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=gpu
#SBATCH --gres gpu:1
#SBATCH --qos long
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com  


method="scglue"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments
parse_arguments "$@"

# Pass arguments to Python script
python_args="--rna $rna --prediction $prediction"
if [ ! -z "$atac" ]; then
    python_args="$python_args --atac $atac"
fi

singularity run --nv ../../images/${method} python src/methods/${method}/script.py $python_args
#!/bin/bash
#SBATCH --job-name=scenic
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=40:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

method="scenic"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments (temp_dir handled separately)
temp_dir=""
filtered_args=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --temp_dir) temp_dir="$2"; shift 2 ;;
        *) filtered_args+=("$1"); shift ;;
    esac
done
parse_arguments "${filtered_args[@]}"

# Derive dataset-specific temp dir to avoid collisions between concurrent runs
if [ -z "$temp_dir" ]; then
    dataset=$(basename $(dirname $prediction))
    temp_dir="output/${method}/${dataset}"
fi

# Pass arguments to Python script
python_args="--prediction $prediction --temp_dir $temp_dir --tf_all ${tf_all:-resources/grn_benchmark/prior/tf_all.csv}"
if [ ! -z "$rna" ]; then
    python_args="$python_args --rna $rna"
fi
if [ ! -z "$layer" ]; then
    python_args="$python_args --layer $layer"
fi

singularity run resources/singularity/scenic python src/methods/${method}/script.py $python_args
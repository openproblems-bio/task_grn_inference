#!/bin/bash
#SBATCH --job-name=positive_control
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

method="positive_control"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments
parse_arguments "$@"

# Pass arguments to Python script
python_args="--prediction $prediction --tf_all ${tf_all:-resources/grn_benchmark/prior/tf_all.csv}"

if [ ! -z "$rna_all" ]; then
    python_args="$python_args --rna_all $rna_all"
fi

if [ ! -z "$layer" ]; then
    python_args="$python_args --layer $layer"
fi

echo $python_args

# Initialize conda in the (non-interactive) batch shell, then activate the env.
# py10 is itself a conda install, so on PATH `conda` resolves to py10's conda whose
# base lacks genernbi. $CONDA_EXE (exported on activation) points to the outer conda;
# use it to source the right conda.sh so `conda activate genernbi` resolves correctly.
source "$("$CONDA_EXE" info --base)/etc/profile.d/conda.sh"
conda activate genernbi
python src/methods/${method}/script.py $python_args
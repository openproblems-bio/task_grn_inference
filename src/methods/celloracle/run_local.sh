#!/bin/bash
#SBATCH --job-name=celloracle
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

NEW_CACHE=$TMPDIR/cache
mkdir -p $NEW_CACHE
if [ -z $XDG_CACHE_HOME ]; then
    XDG_CACHE_HOME=$HOME/.cache
fi
cp -r $XDG_CACHE_HOME/gimmemotifs $NEW_CACHE/
export XDG_CACHE_HOME=$NEW_CACHE
echo 'Using for cache' $XDG_CACHE_HOME 

method="celloracle"

# Import argument parsing functionality
source "src/utils/parse_args.sh"

# Parse command line arguments
parse_arguments "$@"

# Pass arguments to Python script
python_args="--rna $rna --atac $atac --prediction $prediction"
if [ ! -z "$layer" ]; then
    python_args="$python_args --layer $layer"
fi

singularity run ../../images/celloracle python src/methods/${method}/script.py $python_args
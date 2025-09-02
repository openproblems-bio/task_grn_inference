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

dataset=$1
method="celloracle"

if [ -z "$dataset" ]; then
    echo "Error: dataset not provided"
    exit 1
fi

rna="resources/grn_benchmark/inference_data/${dataset}_rna.h5ad"
atac="resources/grn_benchmark/inference_data/${dataset}_atac.h5ad"
prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"

singularity run ../../images/celloracle python src/methods/${method}/script.py --rna $rna --atac $atac --prediction $prediction
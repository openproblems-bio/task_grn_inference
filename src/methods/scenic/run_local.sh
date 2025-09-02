#!/bin/bash
#SBATCH --job-name=scenic
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

dataset=$1
method="scenic"

if [ -z "$dataset" ]; then
    echo "Error: dataset not provided"
    exit 1
fi

rna="resources/grn_benchmark/inference_data/${dataset}_rna.h5ad"
prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"

singularity run ../../images/scenic python src/methods/${method}/script.py --rna $rna --prediction $prediction
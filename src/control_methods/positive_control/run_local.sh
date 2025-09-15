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

dataset=$1
method="positive_control"

if [ -z "$dataset" ]; then
    echo "Error: dataset not provided"
    exit 1
fi

rna_all="resources/extended_data/${dataset}_bulk.h5ad"
prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"

python src/control_methods/${method}/script.py --rna_all $rna_all --prediction $prediction
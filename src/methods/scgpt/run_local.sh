#!/bin/bash
#SBATCH --job-name=scgpt
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=5:00:00
#SBATCH --mem=64GB
#SBATCH --partition=cpu
#SBATCH --qos long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jalil.nourisa@gmail.com

method="scgpt"

source "src/utils/parse_args.sh"
parse_arguments "$@"

singularity run ../../images/${method} python src/methods/${method}/script.py \
    --rna $rna \
    --prediction $prediction \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --temp_dir output/${method}/

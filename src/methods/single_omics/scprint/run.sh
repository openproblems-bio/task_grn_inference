#!/bin/bash
#SBATCH --job-name=scprint
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH --mem=100GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   



# viash run src/methods/single_omics/scprint/config.vsh.yaml -- \
#     --rna resources_test/grn_benchmark/inference_datasets//op_rna.h5ad \
#     --tf_all resources/grn_benchmark/prior/tf_all.csv \
#     --prediction output/prediction.h5ad


# python src/methods/single_omics/scprint/script.py  \
#     --rna resources/grn_benchmark/inference_datasets/op_rna.h5ad \
#     --tf_all resources/grn_benchmark/prior/tf_all.csv \
#     --prediction output/prediction.h5ad

# Exit immediately if a command exits with a non-zero status
set -e
source ~/miniconda3/bin/activate scprint

python src/methods/script_all.py
python src/metrics/script_all_experiment.py
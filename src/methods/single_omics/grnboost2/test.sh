#!/bin/bash
#SBATCH --job-name=grnboos
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

# viash run src/methods/single_omics/grnboost2/config.vsh.yaml -- \
#     --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
#     --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
#     --prediction output/prediction.h5ad \
#     --temp_dir output/grnboost2

singularity run ../../images/scenic python src/methods/single_omics/grnboost2/script.py 
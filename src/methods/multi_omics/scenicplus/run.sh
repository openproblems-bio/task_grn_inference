#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --mem=500GB
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com

# singularity run ../../images/scenicplus python src/methods/multi_omics/scenicplus/script.py 
viash run src/methods/multi_omics/scenicplus/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources/grn_benchmark/inference_data/op_atac.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --prediction output/sp_test/prediction.h5ad \
    --temp_dir output/sp_test/ \
    --num_workers 10
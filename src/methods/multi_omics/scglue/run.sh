#!/bin/bash
#SBATCH --job-name=scglue
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=gpu
#SBATCH --gres gpu:1
#SBATCH --qos long
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com  



viash run src/methods/multi_omics/scglue/config.vsh.yaml -- \
    --rna resources/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources/grn_benchmark/inference_data/op_atac.h5ad \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --prediction output/scglue_new/prediction.h5ad \
    --temp_dir output/scglue_new/ \
    --num_workers 20

# singularity run ../../images/scglue python src/methods/multi_omics/scglue/script.py

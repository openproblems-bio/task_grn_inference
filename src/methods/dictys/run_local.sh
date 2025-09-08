#!/bin/bash
#SBATCH --job-name=dictys
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=gpu
#SBATCH --gres gpu:1
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

images_dir="/home/jnourisa/projs/images"

# singularity run ../../images/dictys \
#     python src/methods/dictys/script.py \
#     --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
#     --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad 

singularity run $images_dir/dictys bash dictys_helper split_bam.sh bams.bam bams --section "CB:Z:" --ref_expression expr.tsv.gz


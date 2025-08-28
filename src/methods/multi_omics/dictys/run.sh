#!/bin/bash
#SBATCH --job-name=dictys
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   


# viash run src/methods/multi_omics/dictys/config.vsh.yaml -- \
#     --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
#     --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
#     --tf_all resources_test/grn_benchmark/prior/tf_all.csv \
#     --prediction output/d_test/prediction.h5ad \
#     --temp_dir output/d_test/ \
#     --num_workers 10



singularity run ../../images/dictys.sif bash src/methods/multi_omics/dictys/tfb.sh \
    --rna_file resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --atac_file resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --output_d output/dictys \
    --input_motif output/d_test/data/motif_file.motif \
    --input_genome resources/extended_data/genome \
    --threads 10


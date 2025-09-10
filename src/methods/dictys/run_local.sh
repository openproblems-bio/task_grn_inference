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

command=(singularity run --nv /home/jnourisa/projs/images/dictys_latest.sif)
# data_dir="output/temp/data/"
# temp_dir="output/dictys_test/"

"${command[@]}" \
    python src/methods/dictys/script.py \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --prediction output/temp/predictions.h5ad \


# "${command[@]}" bash dictys_helper genome_homer.sh hg38 $data_dir/genome

# cd ../makefiles
# "${command[@]}" dictys_helper makefile_template.sh common.mk config.mk env_none.mk static.mk
# singularity run ../../../../../images/dictys.sif bash dictys_helper makefile_update.py ../makefiles/config.mk '{\"DEVICE\": \"cuda:0\", \"GENOME_MACS2\": \"hs\", \"JOINT\": \"1\"}'
# singularity run ../../../../images/dictys.sif bash dictys_helper makefile_check.py

# bash dictys_helper network_inference.sh -j 32 -J 1 static --device cpu

# "${command[@]}" python3 -m dictys  chromatin wellington --nth 4 \
#     $temp_dir/tmp_static/Subset1/reads.bam \
#     $temp_dir/tmp_static/Subset1/reads.bai $temp_dir/tmp_static/Subset1/peaks.bed \
#     $temp_dir/tmp_static/Subset1/footprints.bed
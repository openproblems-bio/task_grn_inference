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

# command=(singularity run --nv /home/jnourisa/projs/images/dictys_latest.sif)
# command=(docker run -it  -v $(pwd):/workspace -w /workspace  ghcr.io/openproblems-bio/task_grn_inference/grn_methods/dictys:dev) 


command=(singularity run --nv /home/jnourisa/projs/external/greta/workflow/envs/dictys.sif)
temp_dir="output/temp/"
data_dir="{$temp_dir}/data/"



# "${command[@]}" python3 -m dictys  chromatin wellington --nth 4 \
#     $temp_dir/tmp_static/cluster_1/reads.bam \
#     $temp_dir/tmp_static/cluster_1/reads.bai $temp_dir/tmp_static/cluster_1/peaks.bed \
#     $temp_dir/tmp_static/cluster_1/footprints.bed

"${command[@]}" \
    python src/methods/dictys/script.py \
    --rna resources_test/grn_benchmark/inference_data/op_rna.h5ad \
    --atac resources_test/grn_benchmark/inference_data/op_atac.h5ad \
    --prediction output/temp/predictions.h5ad \


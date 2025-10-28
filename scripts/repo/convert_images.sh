#!/bin/bash
#SBATCH --job-name=convert_images
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=20:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   


# Control methods
for method in positive_control negative_control pearson_corr; do
    echo "Running control method: $method"
    apptainer exec --no-home --pid \
        docker://ghcr.io/openproblems-bio/task_grn_inference/methods/${method}:build_main \
        /bin/bash -c "echo hi"
done

# GRN inference methods
for method in portia ppcor scenic grnboost2 scenicplus scglue granie figr celloracle scprint; do
    echo "Running GRN method: $method"
    apptainer exec --no-home --pid \
        docker://ghcr.io/openproblems-bio/task_grn_inference/grn_methods/${method}:build_main \
        /bin/bash -c "echo hi"
done
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --job-name=scenicplus
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G 

singularity exec ../../images/scenicplus_old python src/methods/multi_omics/scenicplus/script.py 

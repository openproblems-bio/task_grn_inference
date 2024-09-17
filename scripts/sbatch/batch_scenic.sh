#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=scenic
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G 

singularity exec ../../images/scenic python src/methods/single_omics/scenic/script.py 

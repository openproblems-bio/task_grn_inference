#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=scenic_donor_0
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G 

singularity exec ../../images/scenic python src/methods/single_omics/scenic/script.py 

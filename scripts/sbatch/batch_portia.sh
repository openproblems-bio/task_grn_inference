#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --job-name=portia
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G 

singularity exec ../../images/portia python src/methods/single_omics/portia/script.py 

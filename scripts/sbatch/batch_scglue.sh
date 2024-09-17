#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=scglue_donor_0
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G 

singularity exec ../../images/scglue python src/methods/multi_omics/scglue/script.py 

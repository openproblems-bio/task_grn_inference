#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --mem=64G 
#SBATCH --partition=gpu  
#SBATCH --cpus-per-task=20  
#SBATCH --gres=gpu:1 

singularity exec --nv ../../images/scglue python src/methods/multi_omics/scglue/script.py 

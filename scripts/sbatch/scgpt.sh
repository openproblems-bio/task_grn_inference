#!/bin/bash
#SBATCH --job-name=scgpt
#SBATCH --time=2:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=20  
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# conda activate pytorch
singularity run ../../images/scglue python -c "import torch; print(torch.cuda.is_available())" 

#!/bin/bash
#SBATCH --job-name=skeleton
#SBATCH --time=48:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --mem=128G 
#SBATCH --cpus-per-task=20  
#SBATCH --partition=gpu 
#SBATCH --gres=gpu:1

singularity run ../../images/scglue python src/metrics/skeleton/script.py

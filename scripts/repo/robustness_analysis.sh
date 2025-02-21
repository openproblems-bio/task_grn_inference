#!/bin/bash
#SBATCH --job-name=robustness
#SBATCH --time=48:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --mem=250G 
#SBATCH --cpus-per-task=20 
# SBATCH --partition=gpu 
# SBATCH --gres=gpu:1

python src/robustness_analysis/script_all.py

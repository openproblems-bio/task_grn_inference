#!/bin/bash
#SBATCH --job-name=robustness
#SBATCH --time=48:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=20  

# python src/metrics/script_all.py
python src/robustness_analysis/script_all.py

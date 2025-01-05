#!/bin/bash
#SBATCH --job-name=scores
#SBATCH --time=10:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=1  

# python src/metrics/script_all.py 
python src/metrics/all_metrics/script_all.py 

#!/bin/bash
#SBATCH --job-name=de_analysis
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=500GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   




python src/process_data/de_analysis/script.py 
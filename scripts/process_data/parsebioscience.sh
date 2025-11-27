#!/bin/bash
#SBATCH --job-name=parsebioscience
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=20:00:00
#SBATCH --mem=1500GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   


set -e

echo "Processing bioscience"
python src/process_data/main/parsebioscience/script.py  #--run_test
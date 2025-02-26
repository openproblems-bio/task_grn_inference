#!/bin/bash
#SBATCH --job-name=benchmark_all
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=120GB
#SBATCH --nodelist=bioinf008
#SBATCH --qos long
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

bash scripts/run_benchmark_all.sh
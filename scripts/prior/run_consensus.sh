#!/bin/bash
#SBATCH --job-name=ws
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --time=5-00:00:00
#SBATCH --qos long
#SBATCH --mem=1000GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

set -e

echo "Running consensus for regression 2"
python src/metrics/regression_2/consensus/script.py

echo "Running consensus for ws distance"
python src/metrics/wasserstein/consensus/script.py

echo "Calculating scores for all possible connections, WS distance"
python src/metrics/wasserstein/background_distance/script.py
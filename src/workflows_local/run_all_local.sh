#!/bin/bash
#SBATCH --job-name=run_all_local
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --mem=500GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   


set -e
# - data preprocessing 

# echo "Running opsc_perturbation ..."
# python src/process_data/opsca_perturbation/script.py 

# echo "Running replogle_k562_gwps ..."
# bash src/process_data/replogle_k562_gwps/run.sh 

# echo "Running pereggrn ..."
# bash src/process_data/pereggrn/run.sh 

# echo "Running test_data ..."
# python src/process_data/test_data/datasets/script.py

# - GRN inference and evaluation
echo "Running workflows_local ..."
bash src/workflows_local/benchmark/run.sh


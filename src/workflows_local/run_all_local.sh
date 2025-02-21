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
python src/process_data/opsca_multiomics/script.py 
# python src/process_data/opsca_perturbation/script.py 



# echo "Running pereggrn ..."
# python src/process_data/pereggrn/_script.py
# python src/process_data/replogle/script.py 
# python src/process_data/nakatake/script.py 
# python src/process_data/adamson/script.py 
# python src/process_data/norman/script.py 

# echo "Running replogle_k562_gwps ..." 
# bash src/process_data/replogle_k562_gwps/run.sh # this gets single cell data for replogle

# echo "Running test_data ..."
python src/process_data/test_data/datasets/script.py

aws s3 sync  resources_test/grn_benchmark/ s3://openproblems-data/resources_test/grn/grn_benchmark/ --delete


# - GRN inference and evaluation
# echo "Running workflows_local ..."
# bash src/workflows_local/benchmark/run.sh


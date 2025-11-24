#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=20:00:00
#SBATCH --mem=500GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   


set -e

python src/process_data/main/adamson/script.py 
# python src/process_data/main/nakatake/script.py 
python src/process_data/main/norman/script.py

# echo "Processing opsca"
# python src/process_data/main/opsca/script.py 
# echo "Processing replogle"
# python src/process_data/main/replogle/script.py  #--run_test  #--run_test
# echo "Processing xaira"
# python src/process_data/main/xaira/script.py    #--run_test


# echo "Processing 300BCG"
# python src/process_data/main/300BCG/script.py 
# echo "Processing IBD"
# python src/process_data/main/ibd/script.py 

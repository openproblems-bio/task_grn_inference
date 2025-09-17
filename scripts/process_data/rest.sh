#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   


set -e

# python src/process_data/adamson/script.py 
# python src/process_data/nakatake/script.py 
# python src/process_data/norman/script.py

echo "Processing opsca"
python src/process_data/opsca/script.py 
# echo "Processing replogle"
# python src/process_data/replogle/script.py  #--run_test  #--run_test
# echo "Processing xaira"
# python src/process_data/xaira/script.py    #--run_test
# echo "Processing 300BCG"
# python src/process_data/300BCG/script.py 
# echo "Processing IBD"
# python src/process_data/ibd/script.py 

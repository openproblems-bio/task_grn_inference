#!/bin/bash
#SBATCH --job-name=process_data
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=1000GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   


set -e

# python src/process_data/adamson/script.py 
# python src/process_data/nakatake/script.py 
# python src/process_data/norman/script.py

# python src/process_data/opsca/script.py 
# python src/process_data/replogle/script.py  #--run_test  #--run_test
python src/process_data/xaira/script.py    #--run_test
# python src/process_data/parse_bioscience/script.py  #--run_test

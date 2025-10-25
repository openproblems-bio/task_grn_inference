#!/bin/bash
#SBATCH --job-name=tf_binding_data
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=250GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

python src/process_data/tfb/process_tfb_datasets.py --cell-type hek293
# python src/process_data/tfb/process_tfb_datasets.py --cell-type k562
python src/process_data/tfb/process_tfb_datasets.py --cell-type hct116

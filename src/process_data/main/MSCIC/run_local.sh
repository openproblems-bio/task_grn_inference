#!/bin/bash
#SBATCH --job-name=preprocess_MSCIC
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=128GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jalil.nourisa@gmail.com

singularity exec \
  --bind /vol/projects:/vol/projects \
  resources/singularity/celloracle \
  python3 src/process_data/main/MSCIC/script.py

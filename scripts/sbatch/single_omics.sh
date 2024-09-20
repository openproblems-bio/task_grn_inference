#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G 

echo "singularity exec ../../images/${1} python src/methods/single_omics/${2}/script.py "
singularity exec ../../images/${1} python src/methods/single_omics/${2}/script.py 


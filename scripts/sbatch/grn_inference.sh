#!/bin/bash

#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jalil.nourisa@gmail.com

${1}
# if [ "${1}" = "pearson_corr" ]; then
#     echo "Method: ${1}"
#     python src/control_methods/pearson/script.py ${args}
# elif [ "${1}" = "celloracle" ]; then
#     /home/jnourisa/miniconda3/envs/celloracle/bin/python src/methods/multi_omics/celloracle/script.py ${args}
# elif [ "$1" == "scenic" ] || [ "$1" == "genie3" ] || [ "$1" == "grnboost2" ]; then
#     singularity exec ../../images/scenic python src/methods/single_omics/${1}/script.py ${args}
# elif [ "$1" == "scglue" ]; then
#     singularity exec ../../images/${1} python src/methods/multi_omics/${1}/script.py ${args}
# elif [ "$1" == "ppcor" ]; then #R 
#     singularity exec ../../images/${1} Rscript src/methods/single_omics/${1}/script.R ${args}
# else
#     singularity exec ../../images/${1} python src/methods/single_omics/${1}/script.py ${args}
# fi





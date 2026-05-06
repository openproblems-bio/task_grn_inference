#!/bin/bash
#SBATCH --job-name=scgpt
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=5:00:00
#SBATCH --mem=64GB
#SBATCH --partition=cpu
#SBATCH --qos long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jalil.nourisa@gmail.com

method="scgpt"

source "src/utils/parse_args.sh"
parse_arguments "$@"

# Pre-generate gene annotation cache on host (pybiomart not available in container)
GENE_INFO="resources/grn_benchmark/prior/gene_info_biomart.csv"
if [ ! -f "$GENE_INFO" ]; then
    echo "Generating gene info cache via pybiomart..."
    /home/jnourisa/miniconda3/envs/py10/bin/python -c "
import sys; sys.path.insert(0, 'src')
from utils.util import fetch_gene_info
fetch_gene_info().to_csv('$GENE_INFO')
print('Saved gene info cache to $GENE_INFO')
"
fi

singularity run ../../images/${method} python src/methods/${method}/script.py \
    --rna $rna \
    --prediction $prediction \
    --tf_all resources/grn_benchmark/prior/tf_all.csv \
    --temp_dir output/${method}/

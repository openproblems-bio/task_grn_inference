#!/bin/bash
#SBATCH --job-name=scprint
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

method="scprint"

source "src/utils/parse_args.sh"
parse_arguments "$@"

# Set up a writable lamindb home using the pre-built instance from the container
_lamin_home=$(mktemp -d /tmp/lamin_home_XXXXXX)
_lamin_storage=$(mktemp -d /tmp/lamin_storage_XXXXXX)

singularity exec ../../images/${method}.sif cp -r /root/.lamin/. ${_lamin_home}/.lamin/
singularity exec ../../images/${method}.sif cp /workspace/main/63c3fd677cf055009fc56bed97323c1c.lndb ${_lamin_storage}/

for f in ${_lamin_home}/.lamin/*.env; do
    sed -i "s|lamindb_instance_storage_root=/workspace/main|lamindb_instance_storage_root=${_lamin_storage}|g" "$f"
done

singularity exec --home ${_lamin_home} --bind $(pwd):$(pwd) ../../images/${method}.sif \
    python src/methods/${method}/script.py --rna $rna --prediction $prediction

rm -rf ${_lamin_home} ${_lamin_storage}
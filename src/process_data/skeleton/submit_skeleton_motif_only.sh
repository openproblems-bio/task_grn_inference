#!/bin/bash
# Submit motif-only skeleton building jobs for all non-ATAC datasets.
# Output: resources/grn_benchmark/prior/skeleton_{dataset}.csv
REPO=/home/jnourisa/projs/ongoing/task_grn_inference
SCRIPT=$REPO/src/process_data/skeleton/skeleton_motif_only.py
IMAGE=/home/jnourisa/projs/images/scglue
DATASETS=(nakatake norman replogle parsebioscience 300BCG xaira_HEK293T xaira_HCT116)

for ds in "${DATASETS[@]}"; do
    echo "Submitting skeleton job for: $ds"
    sbatch \
      --job-name="skel_${ds}" \
      --output="$REPO/logs/skel_${ds}_%j.out" \
      --error="$REPO/logs/skel_${ds}_%j.err" \
      --cpus-per-task=8 \
      --mem=32G \
      --time=4:00:00 \
      --partition=cpu \
      --wrap="singularity exec \
        --bind /home/jnourisa/projs:/home/jnourisa/projs \
        $IMAGE \
        python3 $SCRIPT --dataset $ds"
done

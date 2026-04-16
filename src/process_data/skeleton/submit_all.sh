#!/bin/bash
# Submit skeleton building jobs for all datasets.
# Datasets with ATAC: op, ibd_uc, ibd_cd  (full ATAC+motif skeleton, needs GPU + ~48h)
# Datasets without ATAC: use --no_atac (promoter/motif-only, CPU + ~4h)
#
# Usage: bash src/process_data/skeleton/submit_all.sh

REPO=/home/jnourisa/projs/ongoing/task_grn_inference
SCRIPT=$REPO/src/process_data/skeleton/build_skeleton.py
IMAGE=/home/jnourisa/projs/images/scglue
LOGS=$REPO/logs

mkdir -p $LOGS

# ---------------------------------------------------------------------------
# Multiomics datasets (ATAC available) — full ATAC+motif skeleton
# Commented out: skeletons already exist and have been copied to prior/
# Uncomment to rebuild from scratch (requires GPU node, ~48h each).
# ---------------------------------------------------------------------------
#
# DATASETS_ATAC=(op ibd_uc ibd_cd)
#
# for ds in "${DATASETS_ATAC[@]}"; do
#     OUT=$REPO/resources/grn_benchmark/prior/skeleton_${ds}.csv
#     if [ -f "$OUT" ]; then
#         echo "Skipping $ds — already exists: $OUT"
#         continue
#     fi
#     echo "Submitting full ATAC+motif skeleton for $ds"
#     sbatch \
#         --job-name=skel_${ds} \
#         --output=$LOGS/skel_${ds}_%j.out \
#         --error=$LOGS/skel_${ds}_%j.err \
#         --time=48:00:00 \
#         --mem=64G \
#         --cpus-per-task=16 \
#         --gres=gpu:1 \
#         --partition=gpu \
#         --wrap="singularity exec --nv \
#             --bind /home/jnourisa/projs:/home/jnourisa/projs \
#             $IMAGE \
#             python3 $SCRIPT --dataset $ds"
# done

# ---------------------------------------------------------------------------
# Datasets without ATAC — motif-only skeleton
# ---------------------------------------------------------------------------
DATASETS_NO_ATAC=(300BCG nakatake norman parsebioscience replogle xaira_HEK293T xaira_HCT116)

for ds in "${DATASETS_NO_ATAC[@]}"; do
    OUT=$REPO/resources/grn_benchmark/prior/skeleton_${ds}.csv
    if [ -f "$OUT" ]; then
        echo "Skipping $ds — already exists: $OUT"
        continue
    fi
    echo "Submitting motif-only skeleton for $ds"
    sbatch \
        --job-name=skel_${ds} \
        --output=$LOGS/skel_${ds}_%j.out \
        --error=$LOGS/skel_${ds}_%j.err \
        --time=4:00:00 \
        --mem=32G \
        --cpus-per-task=8 \
        --partition=cpu \
        --wrap="singularity exec \
            --bind /home/jnourisa/projs:/home/jnourisa/projs \
            $IMAGE \
            python3 $SCRIPT --dataset $ds --no_atac"
done

echo ""
echo "Existing skeletons (already copied):"
ls $REPO/resources/grn_benchmark/prior/skeleton_*.csv 2>/dev/null

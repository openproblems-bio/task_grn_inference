#!/bin/bash
# Submit skeleton building jobs.
#
# Multiomics datasets (ATAC available): op, ibd_uc, ibd_cd
#   → full ATAC+motif skeleton, dataset-specific, needs GPU + ~48h
#
# Non-ATAC datasets: ONE shared skeleton per gene-ID type (not per dataset)
#   → skeleton_motif.csv         (gene_name: 300BCG, nakatake, norman,
#                                  parsebioscience, replogle)
#   → skeleton_motif_ensembl.csv (gene_id:  xaira_HEK293T, xaira_HCT116)
#
# Usage: bash src/process_data/skeleton/submit_all.sh

REPO=/home/jnourisa/projs/ongoing/task_grn_inference
SCRIPT=$REPO/src/process_data/skeleton/skeleton_motif_only.py
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
# Shared motif-only skeleton — ONE file for ALL non-ATAC datasets
#
#   skeleton_motif.csv — built from tss_h38.bed (all protein-coding genes,
#                        ~19k unique genes). Dataset-agnostic.
# ---------------------------------------------------------------------------

OUT_MOTIF=$REPO/resources/grn_benchmark/prior/skeleton_motif.csv
if [ -f "$OUT_MOTIF" ]; then
    echo "Skipping skeleton_motif.csv — already exists"
else
    echo "Submitting shared motif skeleton (protein-coding genes from tss_h38.bed)"
    sbatch \
        --job-name=skel_motif \
        --output=$LOGS/skel_motif_%j.out \
        --error=$LOGS/skel_motif_%j.err \
        --time=4:00:00 \
        --mem=32G \
        --cpus-per-task=8 \
        --partition=cpu \
        --wrap="singularity exec \
            --bind /home/jnourisa/projs:/home/jnourisa/projs \
            $IMAGE \
            python3 $SCRIPT"
fi

echo ""
echo "Existing skeletons (already copied):"
ls $REPO/resources/grn_benchmark/prior/skeleton_*.csv 2>/dev/null

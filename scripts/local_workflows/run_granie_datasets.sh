#!/bin/bash
# Run GraNIE GRN inference locally on all multiomics datasets (op, ibd_cd, ibd_uc).
# Run from the repo root: bash scripts/local_workflows/run_granie_datasets.sh
#
# Optional args:
#   --datasets   comma-separated dataset list  (default: op,ibd_cd,ibd_uc)
#   --layer      expression layer              (default: lognorm)
#   --dry_run    print commands without running

set -e

DATASETS_ARG="op,ibd_cd,ibd_uc"
LAYER="lognorm"
DRY_RUN=false

for arg in "$@"; do
    case $arg in
        --datasets=*)  DATASETS_ARG="${arg#*=}" ;;
        --layer=*)     LAYER="${arg#*=}" ;;
        --dry_run)     DRY_RUN=true ;;
        *)             echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

IFS=',' read -ra DATASETS <<< "$DATASETS_ARG"

RESOURCES_DIR="resources"
PREDICTIONS_DIR="${RESOURCES_DIR}/results"
METHOD="granie"

echo "=========================================="
echo "GraNIE local inference"
echo "Datasets : ${DATASETS[*]}"
echo "Layer    : ${LAYER}"
echo "Dry run  : ${DRY_RUN}"
echo "=========================================="

for dataset in "${DATASETS[@]}"; do
    rna="${RESOURCES_DIR}/grn_benchmark/inference_data/${dataset}_rna.h5ad"
    atac="${RESOURCES_DIR}/grn_benchmark/inference_data/${dataset}_atac.h5ad"
    rna_all="${RESOURCES_DIR}/extended_data/${dataset}_rna_all.h5ad"
    out_dir="${PREDICTIONS_DIR}/${dataset}"
    prediction="${out_dir}/${dataset}.${METHOD}.${METHOD}.prediction.h5ad"

    echo ""
    echo "--- Dataset: ${dataset} ---"

    # Validate inputs
    for f in "$rna" "$atac"; do
        if [[ ! -f "$f" ]]; then
            echo "  [SKIP] Missing input: $f"
            continue 2
        fi
    done

    mkdir -p "$out_dir"

    cmd="bash src/methods/${METHOD}/run_local.sh \
        --rna $rna \
        --atac $atac \
        --rna_all $rna_all \
        --prediction $prediction \
        --layer $LAYER"

    echo "  Output: $prediction"
    echo "  CMD: $cmd"

    if [[ "$DRY_RUN" == "false" ]]; then
        eval "$cmd"
        echo "  [DONE] $dataset"
    else
        echo "  [DRY RUN] skipping execution"
    fi
done

echo ""
echo "All done."

#!/bin/bash
#SBATCH --job-name=ws
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --time=10:00:00
#SBATCH --qos long
#SBATCH --mem=1000GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

set -e

DATASET=""
NEW_MODEL_PATH=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --dataset)
      DATASET="$2"
      shift 2
      ;;
    --new_model)
      NEW_MODEL_PATH="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: sbatch run_consensus.sh --dataset <dataset> [--new_model <path>]"
      exit 1
      ;;
  esac
done

if [ -z "$DATASET" ]; then
  echo "Usage: sbatch run_consensus.sh --dataset <dataset> [--new_model <path>]"
  exit 1
fi

models_dir="resources/results/$DATASET"
models=("pearson_corr" "positive_control" "portia" "ppcor" "scenic" "scprint" "grnboost" "scenicplus" "scglue" "granie" "figr" "celloracle" "scgpt" "geneformer" "spearman_corr")
python src/utils/config.py
source src/utils/config.env
METHODS=(${METHODS//,/ })

predictions=()
for model in "${METHODS[@]}"; do
    file="${models_dir}/${DATASET}.${model}.${model}.prediction.h5ad"
    if [ -e "$file" ]; then
        predictions+=("$file")
    fi
done

if [ -n "$NEW_MODEL_PATH" ]; then
    if [ -e "$NEW_MODEL_PATH" ]; then
        echo "Adding new model: $NEW_MODEL_PATH"
        predictions+=("$NEW_MODEL_PATH")
    else
        echo "Warning: New model path does not exist: $NEW_MODEL_PATH"
    fi
fi

printf '%s\n' "${predictions[@]}"

echo "Running consensus for Regression"
datasets=(${DATASET})
for dataset in "${datasets[@]}"; do
    echo "Running regression consensus for dataset: $dataset"
    
    python src/metrics/regression/consensus/script.py \
        --dataset "$dataset" \
        --regulators_consensus "resources/grn_benchmark/prior/regulators_consensus_${dataset}.json" \
        --evaluation_data "resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad" \
        --predictions "${predictions[@]}"
done


echo "Running consensus for ws distance"
applicable_datasets=("norman" "adamson" "replogle" "xaira_HEK293T" "xaira_HCT116" )
for dataset in "${datasets[@]}"; do
    skip=true
    for d in "${applicable_datasets[@]}"; do
        if [[ "$dataset" == "$d" ]]; then
            skip=false
            break
        fi
    done

    if $skip; then
        echo "Skipping dataset: $dataset"
        continue
    fi

    echo "Running ws distance consensus for dataset: $dataset"
    python src/metrics/ws_distance/consensus/script.py \
        --dataset "$dataset" \
        --models_dir "$models_dir" \
        --ws_consensus "resources/grn_benchmark/prior/ws_consensus_${dataset}.csv" \
        --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
        --evaluation_data_sc "resources/processed_data/${dataset}_evaluation_sc.h5ad" \
        --models "${models[@]}"
done


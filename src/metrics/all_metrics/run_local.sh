#!/bin/bash
#SBATCH --job-name=metrics
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2:00:00
#SBATCH --mem=120GB
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

set -euo pipefail

layer="lognorm" 
reg_type="ridge"
num_workers=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset)
      dataset="$2"
      shift 2
      ;;
    --prediction)
      prediction="$2"
      shift 2
      ;;
    --layer)
      layer="$2"
      shift 2
      ;;
    --score)
      score="$2"
      shift 2
      ;;
    --reg_type)
      reg_type="$2"
      shift 2
      ;;
    --num_workers)
      num_workers="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done
if [[ -n "$layer" ]]; then
  echo "Layer is set to: $layer"
fi

# Check required args
if [[ -z "${dataset:-}" || -z "${prediction:-}" || -z "${score:-}" ]]; then
  echo "Usage: $0 --dataset <name> --prediction <file> --score <metric>"
  exit 1
fi


# Run metrics
python src/metrics/all_metrics/script.py \
  --prediction "${prediction}" \
  --evaluation_data "resources/grn_benchmark/evaluation_data/${dataset}_bulk.h5ad" \
  --evaluation_data_sc "resources/grn_benchmark/evaluation_data/${dataset}_sc.h5ad" \
  --evaluation_data_de "resources/grn_benchmark/evaluation_data/${dataset}_de.h5ad" \
  --regulators_consensus "resources/grn_benchmark/prior/regulators_consensus_${dataset}.json" \
  --layer "${layer}" \
  --reg_type "${reg_type}" \
  --num_workers "${num_workers}" \
  --tf_all "resources/grn_benchmark/prior/tf_all.csv" \
  --ws_consensus "resources/grn_benchmark/prior/ws_consensus_${dataset}.csv" \
  --ws_distance_background "resources/grn_benchmark/prior/ws_distance_background_${dataset}.csv" \
  --score "${score}"
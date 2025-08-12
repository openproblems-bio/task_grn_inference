set -e


bash scripts/run_grn_inference.sh
# bash scripts/prior/run_consensus.sh # run consensus for regression 2 and ws distance -> needs to be run after adding each method and dataset
# bash scripts/prior/run_ws_background.sh # run background distance for ws distance -> needs to be run after adding each dataset
bash scripts/run_grn_evaluation.sh
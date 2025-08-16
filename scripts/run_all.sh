set -e

datasets=('parsebioscience' 'xaira_HCT116' 'xaira_HEK293T') #'replogle' 'op' 'nakatake' 'adamson' 'norman' 
run_local=false

run_grn_inference=false
run_grn_evaluation=false
run_download=true

for dataset in "${datasets[@]}"; do
    if [ "$run_grn_inference" = true ]; then
        echo "Running GRN inference for dataset: $dataset"
        if [ "$run_local" = true ]; then
            echo "Running locally"
        else
            echo "Running on AWS"
        fi
        bash scripts/run_grn_inference.sh $dataset $run_local
        
    fi

    if [ "$run_grn_evaluation" = true ]; then
        if [ "$run_local" = false ]; then
            echo "Downloading inference results from AWS"
            aws s3 sync  s3://openproblems-data/resources/grn/results/$dataset resources/results/$dataset 
        fi 
        echo "Running consensus for dataset: $dataset"
        bash scripts/prior/run_consensus.sh $dataset # run consensus for regression 2 and ws distance -> needs to be run after adding each method and dataset
        
        if [ "$run_local" = false ]; then
            echo "Syncing prior results to AWS"
            aws s3 sync  resources/grn_benchmark/prior s3://openproblems-data/resources/grn/prior 
        fi

        echo "Running GRN evaluation for dataset: $dataset"
        bash scripts/run_grn_evaluation.sh $dataset $run_local
    fi

    if [ "$run_download" = true ]; then
        if [ "$run_local" = false ]; then
            echo "Downloading evaluation results from AWS"
            aws s3 sync  s3://openproblems-data/resources/grn/results/$dataset resources/results/$dataset 
        fi
    fi
    # bash scripts/prior/run_ws_background.sh # run background distance for ws distance -> needs to be run after adding each dataset

done





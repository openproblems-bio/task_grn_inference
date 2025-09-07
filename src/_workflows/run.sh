
#!/bin/bash
set -e

run_prefix='sbatch' #bash
DATASETS=('xaira_HCT116' ) #'op' 'adamson' 'replogle' 'norman' 'nakatake' 'parsebioscience'  '300BCG' 'xaira_HCT116' 'xaira_HEK293T'
# METHODS=('negative_control' 'positive_control' 'pearson_corr' 'portia' 'ppcor' 'grnboost' 'scenic'  'scenicplus' 'scglue' 'figr' 'granie')
METHODS=( 'grnboost' 'scenic')

methods_dir='src/methods/'
ctr_methods_dir='src/control_methods/'



run_func() {
    local method=$1
    local dataset=$2

    echo "Running $method on $dataset"

    local script
    if [[ "$method" == "negative_control" || "$method" == "positive_control" || "$method" == "pearson_corr" ]]; then
        script="${ctr_methods_dir}${method}/run_local.sh"
    else
        script="${methods_dir}${method}/run_local.sh"
    fi

    if [[ "$run_prefix" == "bash" ]]; then
        bash "$script" "$dataset"
    elif [[ "$run_prefix" == "sbatch" ]]; then
        # submit the job and capture the job ID
        output=$(sbatch "$script" "$dataset")
        echo "$output"
        # sbatch usually returns: "Submitted batch job 12345678"
        jobid=$(echo "$output" | awk '{print $4}')
        echo "Job ID: $jobid"
    else
        echo "Unknown run_prefix: $run_prefix"
    fi
}

for dataset in "${DATASETS[@]}"; do
    mkdir resources/results/$dataset || true
    for method in "${METHODS[@]}"; do
        run_func "$method" "$dataset"
    done
done
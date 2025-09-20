
#!/bin/bash
set -e

run_prefix='bash' #bash
DATASETS=('xaira_HEK293T' ) #'op' 'adamson' 'replogle' 'norman' 'nakatake' 'parsebioscience'  '300BCG' 'xaira_HCT116' 'xaira_HEK293T'
# METHODS=('negative_control' 'positive_control' 'pearson_corr' 'portia' 'ppcor' 'grnboost' 'scenic'  'scenicplus' 'scglue' 'figr' 'granie')
METHODS=( 'positive_control' 'pearson_corr')

methods_dir='src/methods/'
ctr_methods_dir='src/methods/'



run_func() {
    local method=$1
    local dataset=$2

    echo "Running $method on $dataset"

    script="${methods_dir}${method}/run_local.sh"

    rna="resources/grn_benchmark/inference_data/${dataset}_rna.h5ad"
    rna_all="resources/extended_data/${dataset}_bulk.h5ad"
    atac="resources/grn_benchmark/inference_data/${dataset}_atac.h5ad"
    prediction="resources/results/${dataset}/${dataset}.${method}.${method}.prediction.h5ad"

    arguments="--rna $rna --prediction $prediction --atac $atac --rna_all $rna_all"

    if [[ "$run_prefix" == "bash" ]]; then
        bash "$script" $arguments
    elif [[ "$run_prefix" == "sbatch" ]]; then
        # submit the job and capture the job ID
        output=$(sbatch "$script" $arguments)
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
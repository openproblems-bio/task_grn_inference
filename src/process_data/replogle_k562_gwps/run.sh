#!/bin/bash
#SBATCH --job-name=run_replogle
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00
#SBATCH --mem=500GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

# python src/process_data/replogle_k562_gwps/sparsify_raw_counts.py  #one time at the begining to covert the downloaded data to the sparsified raw counts

python src/process_data/replogle_k562_gwps/script.py \
        --input resources/datasets_raw/replogle_K562_gwps_raw_singlecell.h5ad \
        --tf_all resources/grn_benchmark/prior/tf_all.csv  \
        --adata_test_sc resources/grn_benchmark/evaluation_data/replogle_sc.h5ad \
        --adata_train_sc resources/extended_data/replogle_train_sc.h5ad \
        --adata_train_sc_subset resources/grn_benchmark/inference_data/replogle_rna_sc_subset.h5ad
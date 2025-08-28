#!/bin/bash
#SBATCH --job-name=sync_resources
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL      
#SBATCH --mail-user=jalil.nourisa@gmail.com   

set -e

# aws s3 sync  s3://openproblems-data/resources/grn/grn_benchmark resources/grn_benchmark  --no-sign-request


# aws s3 sync  s3://openproblems-data/resources/grn/grn_models resources/grn_models --delete 
# aws s3 sync  resources_test/ s3://openproblems-data/resources_test/grn/  --delete 


# aws s3 sync  s3://openproblems-data/resources/grn/extended_data   resources/extended_data/

# aws s3 sync  resources/grn_benchmark/ s3://openproblems-data/resources/grn/grn_benchmark  
# aws s3 sync  resources/extended_data/ s3://openproblems-data/resources/grn/extended_data  
# aws s3 sync resources/results/  s3://openproblems-data/resources/grn/results/ --delete 
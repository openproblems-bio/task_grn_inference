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

# Push results to S3, skipping files larger than 1GB
EXCLUDES=$(find resources/results/ -type f -size +1G \
  | sed 's|resources/results/||' \
  | awk '{print "--exclude \"" $0 "\""}' \
  | tr '\n' ' ')
eval aws s3 sync resources/results/ s3://openproblems-data/resources/grn/results/ $EXCLUDES --delete

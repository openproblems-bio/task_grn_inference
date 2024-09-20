#!/bin/bash

# Define job IDs and corresponding job names
job_ids=("7742241" "7742249" "7742283" "7742285" "7742364" "7742343")
job_names=("portia" "grnboost2" "scenic" "genie3" "ppcor" "scglue")

# Loop through each job ID and job name simultaneously
for i in "${!job_ids[@]}"
do
    job_id=${job_ids[$i]}
    job_name=${job_names[$i]}
    
    echo "Fetching resource usage for Job: $job_name (ID: $job_id)"
    
    # Run the sacct command and write the output to a file using the job name instead of job ID
    sacct -j "$job_id" --format=JobID,JobName,AllocCPUs,Elapsed,State,MaxRSS,MaxVMSize > "resources/results/d0_hvgs/${job_name}.txt"
    
    echo "Resource usage for Job $job_name saved to resources/results/d0_hvgs/${job_name}.txt"
done
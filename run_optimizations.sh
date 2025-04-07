#!/bin/bash

# Load Conda module (if required)
module load conda

# Activate your Conda environment
conda activate myenv

# Get the total number of parameter combinations from R
NUM_COMBOS=$(Rscript -e 'Lambda <- readRDS("opt_results/data/Lambda.rds"); B <- readRDS("opt_results/data/B.rds"); cat(nrow(expand.grid(lambda=Lambda, beta=B)))')

echo "Starting optimization for $NUM_COMBOS parameter combinations..."

# Create results directory if it doesn't exist
mkdir -p opt_results/oracular/indiv_res

# Loop over each parameter index and submit a job for each optimization
for ((i=1; i<=NUM_COMBOS; i++)); do
  echo "[$(date +'%H:%M:%S')] Submitting job for optimization index $i..."
  
  # Submit the job using oarsub (requesting 12 cores per job)
  oarsub -l "host=1/core=12" -n "opt_job_$i" "Rscript run_optimization.R $i"
done

echo "All optimization jobs submitted."

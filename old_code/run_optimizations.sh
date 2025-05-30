#!/bin/bash

# Get the total number of parameter combinations from R
NUM_COMBOS=$(Rscript -e 'Lambda <- readRDS("opt_results/data/Lambda.rds"); B <- readRDS("opt_results/data/B.rds"); cat(nrow(expand.grid(lambda=Lambda, beta=B)))')

echo "Starting optimization for $NUM_COMBOS parameter combinations..."

# Create results directory if it doesn't exist
mkdir -p opt_results/oracular/indiv_res

# Loop over each parameter index and submit a job for each optimization
for ((i=1; i<=NUM_COMBOS; i++)); do
  echo "[$(date +'%H:%M:%S')] Submitting job for optimization index $i..."
  # Submit the job using oarsub (requesting 12 cores per job)
  oarsub -l "host=1/core=12"\
    -O /dev/null -E /dev/null \
    "module load conda && conda activate myenv && Rscript run_optimization.R $i"
done
oarsub -l "host=1/core=12" \
  -O /dev/null -E /dev/null\
  "module load conda && conda activate myenv && Rscript load_results.R"

echo "All optimization jobs submitted."

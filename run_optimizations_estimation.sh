#!/bin/bash

# Get the total number of parameter combinations from R
NUM_COMBOS=$(Rscript -e 'Lambda <- readRDS("opt_results/data/Lambda.rds"); B <- readRDS("opt_results/data/B.rds"); cat(nrow(expand.grid(lambda=Lambda, beta=B)))')

echo "Starting optimization for $NUM_COMBOS parameter combinations..."

# Create results directory if it doesn't exist
mkdir -p opt_results/oracular/indiv_res

# Loop over each parameter index and submit a job for each optimization
for ((i=1; i<=NUM_COMBOS; i++)); do
  echo "[$(date +'%H:%M:%S')] Submitting job for optimization index $i..."
  OUT_LOG="opt_results/logs/opt_job_${i}_out.stdout"
  ERR_LOG="opt_results/logs/opt_job_${i}_err.stderr"
  # Submit the job using oarsub (requesting 12 cores per job)
  oarsub -l "host=1/core=12" -n "opt_job_$i" \
    "module load conda && conda activate myenv && Rscript run_optimization.R $i > $OUT_LOG 2> $ERR_LOG"
done
OUT_LOG="opt_results/logs/Final_output_out.stdout"
ERR_LOG="opt_results/logs/Final_output_err.stderr"
oarsub -l "host=1/core=12" -n  "Final_output"\
 "module load conda && conda activate myenv && Rscript load_results.R $i > $OUT_LOG 2> $ERR_LOG"

echo "All optimization jobs submitted."

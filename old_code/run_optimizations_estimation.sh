#!/bin/bash

# Max number of jobs allowed to wait in queue
BATCH_SIZE=100

# Function to get current number of user's pending jobs
get_job_count() {
  oarstat -u $USER | grep -c "Waiting"
}

# Get total number of parameter combinations
NUM_COMBOS=$(Rscript -e '
  param_grid <- readRDS("results/data/grid_est.rds"); 
  cat(nrow(param_grid))')

for ((i=1; i<=NUM_COMBOS; i++)); do
  echo "[$(date +'%H:%M:%S')] Submitting job for optimization index $i..."
  
  # Submit job
  oarsub -l "host=1/core=12"\
  -O /dev/null -E /dev/null\
  "module load conda && conda activate myenv && Rscript run_optimization_est.R $i"
  
  # Wait for the batch to finish before submitting the next batch
  if (( i % BATCH_SIZE == 0 )); then
    echo "[$(date +'%H:%M:%S')] Sleeping for a while to let jobs clear from the queue."
    sleep 400  # Wait for 5 minutes (adjust if necessary)
  fi
done

# Final job to combine/load results
oarsub -l "host=1/core=12" \
  -O /dev/null -E /dev/null \
  "module load conda && conda activate myenv && Rscript load_results_est.R"

echo "All optimization jobs submitted."

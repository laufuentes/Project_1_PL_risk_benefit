#!/bin/bash

# Get the total number of parameter combinations from R
B=$(Rscript -e 'cat(readRDS(file.path("MC","results","data","MC_iter.rds")))')

echo "Starting nuissance parameter computation for $B datasets..."

# Loop over each parameter index and submit a job for each optimization
for ((i=1; i<=B; i++)); do
  echo "Submitting job for nuisance parameter estimation index $i..."
  oarsub -l "host=1/core=8"\
    -O /dev/null -E /dev/null \
    "module load conda && conda activate myenv && Rscript MC/compute_nuis_params.R $i"
done
echo "All estimation jobs submitted."
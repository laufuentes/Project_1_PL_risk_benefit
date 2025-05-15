#!/bin/bash

folder="Simple"

# while read -r line; do
#   mkdir -p "$folder/$line"
# done < run_MC_simu/folder_structure.txt

# # Get the number of Monte Carlo iterations (B)
# echo "Starting data generation process ..."

# # Submit the job to generate scenarios
# oarsub -l "host=1/core=8" \
#     -O /dev/null -E /dev/null \
#     "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/generate_scenarios.R $folder"

# # Wait for the data generation to (hopefully) finish
# echo "Waiting for data generation to complete..."
# sleep 300

B=$(Rscript -e "cat(readRDS(file.path('$folder', 'results', 'data', 'MC_iter.rds')))")
echo "Starting nuisance parameter computation for $B datasets in scenario '$folder'..."
# Loop over each MC iteration and submit jobs for nuisance parameter estimation
for ((i=1; i<=B; i++)); do
  echo "Submitting job for nuisance parameter estimation index $i..."
  oarsub -l "host=1/core=8" \
      -O /dev/null -E /dev/null \
      "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/compute_nuis_params.R $i"
done

echo "All estimation jobs submitted."
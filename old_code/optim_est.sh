#!/bin/bash

# Get counts from R
B=$(Rscript -e 'cat(readRDS(file.path("MC", "results", "data", "MC_iter.rds")))')
L=$(Rscript -e 'cat(nrow(readRDS(file.path("MC", "results", "data", "param_combinations.rds"))))')

echo "Starting nuisance parameter computation for $B datasets and $L parameter combinations..."

max_jobs=199

for ((i=1; i<=B; i++)); do
  for ((l=1; l<=L; l++)); do
    submitted=0
    while [ $submitted -eq 0 ]; do
      current_jobs=$(oarstat -u $USER | grep -cE "Running|Waiting")

      if (( current_jobs < max_jobs )); then
        output=$(oarsub -l "host=1/core=12" \
          -O /dev/null -E /dev/null \
          "module load conda && conda activate myenv && Rscript MC/optim_comb_multiple_folds.R $i $l")

        if echo "$output" | grep -q "OAR_JOB_ID"; then
          echo "Job submitted for i=$i, l=$l → $output"
          submitted=1
        else
          echo "oarsub failed to submit (possibly due to limit). Waiting 6 mins..."
          sleep 360
        fi
      else
        echo "[$(date +'%H:%M:%S')] Queue full ($current_jobs jobs). Waiting 6 mins..."
        sleep 360
      fi
    done
  done
done

echo "All estimation jobs submitted."

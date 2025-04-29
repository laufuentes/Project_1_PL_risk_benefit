#!/bin/bash

# Get counts from R
folder=$(Rscript -e 'cat(readRDS(file.path("run_MC_simu", "current_scenario.rds")))')
# Now get B from inside that folder
B=$(Rscript -e "cat(readRDS(file.path('${folder}', 'results', 'data', 'MC_iter.rds')))")
L=$(Rscript -e "cat(nrow(readRDS(file.path('${folder}', 'results', 'data', 'param_combinations.rds'))))")

echo "Starting nuisance parameter computation for $B datasets and $L parameter combinations..."
max_jobs=100

# for ((i=1; i<=B; i++)); do
#   for ((l=1; l<=L; l++)); do
#     submitted=0
#     while [ $submitted -eq 0 ]; do
#       current_jobs=$(oarstat -u $USER | awk 'NR>1' | wc -l)

#       if (( current_jobs < max_jobs )); then
#         output=$(oarsub -q "abaca" -l "core=8" \
#           -O /dev/null -E /dev/null \
#           "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/optim_comb.R $i $l")

#         if echo "$output" | grep -q "OAR_JOB_ID"; then
#           echo "Job submitted for i=$i, l=$l → $output"
#           submitted=1
#         else
#           echo "oarsub failed (maybe too busy?). Waiting 3 min..."
#           sleep 180
#         fi
#       else
#         echo "[$(date +'%H:%M:%S')] Too many jobs ($current_jobs). Waiting 3 min..."
#         sleep 180
#       fi
#     done
#   done
# done

# echo "All optimization jobs submitted."
# sleep 180
echo "Starting optimization-estimation for $B datasets and $L parameter combinations..."



for ((i=1; i<=B; i++)); do
  for ((l=1; l<=L; l++)); do
    submitted=0
    while [ $submitted -eq 0 ]; do
      current_jobs=$(oarstat -u $USER | awk 'NR>1' | wc -l)

      if (( current_jobs < max_jobs )); then
        output=$(oarsub -q "abaca" -l "core=8" \
          -O /dev/null -E /dev/null \
          "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/optim_comb_multiple_folds.R $i $l")

        if echo "$output" | grep -q "OAR_JOB_ID"; then
          echo "Job submitted for i=$i, l=$l → $output"
          submitted=1
        else
          echo "oarsub failed (maybe too busy?). Waiting 3 min..."
          sleep 180
        fi
      else
        echo "[$(date +'%H:%M:%S')] Too many jobs ($current_jobs). Waiting 3 min..."
        sleep 180
      fi
    done
  done
done

echo "All estimation jobs submitted."

i=1
l=1
j=2
K=200 
folder=$(Rscript -e 'cat(readRDS(file.path("run_MC_simu", "current_scenario.rds")))')
# Function: Extract job ID from oarsub output
get_job_id_from_output() {
    echo "$1" | grep -oP 'OAR_JOB_ID=\K\d+'
}

# Function: Wait for job to finish based on actual status
wait_for_job_completion() {
    local job_id="$1"
    echo "Waiting for job $job_id to finish..."

    while true; do
        job_info=$(oarstat -f -j "$job_id" 2>/dev/null)

        if [[ -z "$job_info" ]]; then
            echo "Job $job_id not visible yet. Sleeping 10s..."
            sleep 10
            continue
        fi

        job_status=$(echo "$job_info" | grep -i "state" | head -1 | awk '{print $NF}')
        echo "Job $job_id status: $job_status"

        if [[ "$job_status" == "Terminated" || "$job_status" == "Error" || "$job_status" == "Canceled" ]]; then
            echo "Job $job_id finished with status: $job_status"
            break
        fi

        sleep 30
    done
}

echo "Starting sequential job submission: i=$i, l=$l, j=$j"

# ---------- INITIAL JOB ----------
echo "Submitting initial job..."
output=$(oarsub -q abaca -l core=12,walltime=01:00:00 \
    -O /dev/null -E /dev/null \
    "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/alternation-procedure/optim_init.R $i $l $j")

echo "OARSUB OUTPUT: $output"

job_id=$(get_job_id_from_output "$output")
if [[ -z "$job_id" ]]; then
    echo "ERROR: Failed to extract job ID from initial job submission"
    exit 1
fi
echo "Initial job submitted with job ID: $job_id"

# ---------- WAIT ----------
wait_for_job_completion "$job_id"

# ---------- ITERATIVE JOBS ----------
for ((k=1;k<=K;k++)); do
    echo "-----------------------------------------------------"

    # Check for final result file
    final_file="${folder}/results/new_estimated/${i}${l}${j}_final_theta.rds"
    if [[ -f "$final_file" ]]; then
        echo "Final result file '$final_file' found. Stopping further submissions."
        break
    fi

    echo "Submitting job for iteration $k..."
    output=$(oarsub -q abaca -l core=12,walltime=01:00:00 \
        -O /dev/null -E /dev/null \
        "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/alternation-procedure/optim_step.R $i $l $j $k")

    echo "OARSUB OUTPUT: $output"
    
    job_id=$(get_job_id_from_output "$output")
    if [[ -z "$job_id" ]]; then
        echo "ERROR: Failed to extract job ID for iteration $k"
        exit 1
    fi
    echo "Job $k submitted with ID: $job_id"

    # ---------- WAIT ----------
    wait_for_job_completion "$job_id"
done

echo "All jobs completed sequentially or were stopped due to final file presence."




# folder=$(Rscript -e 'cat(readRDS(file.path("run_MC_simu", "current_scenario.rds")))')
# # Now get B from inside that folder
# B=$(Rscript -e "cat(readRDS(file.path('${folder}', 'results', 'data', 'MC_iter.rds')))")
# L=$(Rscript -e "cat(nrow(readRDS(file.path('${folder}', 'results', 'data', 'param_combinations.rds'))))")
# J=$(Rscript -e "cat(readRDS(file.path('${folder}', 'results', 'data', 'Jfold.rds')))")

# echo "folder: $folder"
# echo "B: $B"
# echo "L: $L"
# echo "J: $J"

# echo "Starting nuisance parameter computation for $B datasets and $L parameter combinations..."
# max_jobs=100

# for ((i=1; i<=1; i++)); do
#     for ((l=1; l<=L; l++)); do
#         for ((j=1; j<=J; j++));do
#             submitted=0
#             while [ $submitted -eq 0 ]; do
#             current_jobs=$(oarstat -u $USER | awk 'NR>1' | wc -l)

#             if (( current_jobs < max_jobs )); then
#                 output=$(oarsub -l "walltime=8:00:00"\
#                  "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/opt_new_est_comb.R $i $l $j")

#                 if echo "$output" | grep -q "OAR_JOB_ID"; then
#                 echo "Job submitted for i=$i, l=$l. j=$j → $output"
#                 submitted=1
#                 else
#                 echo "oarsub failed (maybe too busy?). Waiting 3 min..."
#                 sleep 180
#                 fi
#             else
#                 echo "[$(date +'%H:%M:%S')] Too many jobs ($current_jobs). Waiting 3 min..."
#                 sleep 180
#             fi
#             done
#         done
#     done
# done

#!/bin/bash

# Define initial variables
#!/bin/bash

#!/bin/bash
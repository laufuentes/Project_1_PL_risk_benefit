# Get counts from R
folder=$(Rscript -e 'cat(readRDS(file.path("run_MC_simu", "current_scenario.rds")))')
# Now get B from inside that folder
B=$(Rscript -e "cat(readRDS(file.path('${folder}', 'results', 'data', 'MC_iter.rds')))")
echo "Starting nuisance parameter computation for $B datasets"

oarsub -l "host=1/core=12" \
          -O /dev/null -E /dev/null \
          "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/load_res.R"
sleep 180
for ((i=1; i<=B; i++)); do
    output=$(oarsub -l "host=1/core=12" \
          -O /dev/null -E /dev/null \
          "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/opt_policy_est.R $i")
    echo "Job submitted for i=$i,→ $output"
done 

echo "All estimation jobs submitted."

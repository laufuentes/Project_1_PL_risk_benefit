B=$(Rscript -e 'cat(readRDS(file.path("MC", "results", "data", "MC_iter.rds")))')
echo "Starting nuisance parameter computation for $B datasets"

for ((i=1; i<=B; i++)); do
    output=$(oarsub -l "host=1/core=12" \
          -O /dev/null -E /dev/null \
          "module load conda && conda activate myenv && Rscript MC/opt_policy_est.R $i")
    echo "Job submitted for i=$i,→ $output"
done 

echo "All estimation jobs submitted."

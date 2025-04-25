oarsub -l "host=1/core=12" \
          -O /dev/null -E /dev/null \
          "module load conda && conda activate myenv && Rscript run_MC_simu/R_files/plots.R"
#!/bin/bash
OAR -q production \
    -l host=1/gpu=1 \
    -l walltime=15:00:00 \
    -O OAR_%jobid%.out \
    -E OAR_%jobid%.err 

module load conda
conda activate myenv
Rscript run.R

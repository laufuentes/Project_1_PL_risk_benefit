#!/bin/bash
#oarsub -l host=1/core=12,walltime=2:00:00

module load conda
conda activate myenv
Rscript oracular_simulations.R
#!/bin/bash
#SBATCH --job-name=4lcomp
#SBATCH --output=logs/%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=144:00:00

source /home/xxx/miniconda3/etc/profile.d/conda.sh
conda activate rep2si_env

Rscript /home/xxx/rep2si/comprehension.R

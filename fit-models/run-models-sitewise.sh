#!/bin/bash
#SBATCH --job-name=sitewise
#SBATCH --output=logs/%j.out
#SBATCH --array=1-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=144:00:00

export USE_SUBSET=0
export MODELS_FILE=v3_models_sitewise

source /home/xxx/miniconda3/etc/profile.d/conda.sh
conda activate rep2si_env

srun --label Rscript /home/xxx/rep2si/fit-models-winston-sitewise.R

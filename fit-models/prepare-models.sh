source /home/xxx/miniconda3/etc/profile.d/conda.sh
conda activate rep2si_env

Rscript /home/xxx/rep2si/prepare-model-v2.R
Rscript /home/xxx/rep2si/prepare-model-v3.R
Rscript /home/xxx/rep2si/prepare-model-v3-sitewise.R

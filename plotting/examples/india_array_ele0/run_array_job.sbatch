#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --job-name=india_ele0
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --array=0-9

DATA_DIR=/gws/nopw/j04/ncas_radar_vol2/pestdar/india/raw_data/

DATES=(${DATA_DIR}/*)

DATE_PATH=${DATES[$SLURM_ARRAY_TASK_ID]}


source activate /home/users/eeslb/.bashrc
source activate /home/users/eeslb/miniconda3/envs/pyart_3_8_biodar

python plot_radar_jpr.py ${DATE_PATH} 0 100 /gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/plots/india/

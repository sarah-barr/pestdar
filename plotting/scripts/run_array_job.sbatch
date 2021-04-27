#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --job-name=salalah_test
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=5:00
#SBATCH --array=0-42

DATA_DIR=/gws/nopw/j04/ncas_radar_vol2/pestdar/oman/salalah/raw/cfradial/

DATES=(${DATA_DIR}/*)

DATE_PATH=${DATES[$SLURM_ARRAY_TASK_ID]}

source activate /home/users/eeslb/.bashrc
source activate /home/users/eeslb/miniconda3/envs/pyart_3_8_biodar

python plot_radar_all.py ${DATE_PATH} 1 100 /gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/plots/oman/salalah 

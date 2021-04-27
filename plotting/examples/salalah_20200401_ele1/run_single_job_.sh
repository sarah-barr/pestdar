#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --job-name=salalah_check_mem
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out
#SBATCH --time=00:30:00

source activate /home/users/eeslb/.bashrc
source activate /home/users/eeslb/miniconda3/envs/pyart_3_8_biodar

python plot_radar_all.py /gws/nopw/j04/ncas_radar_vol2/pestdar/oman/salalah/raw/cfradial/20200401/ 1 100 /gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/plots/oman/salalah

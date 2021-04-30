#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --job-name=delhi05
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out
#SBATCH --time=24:00:00

source activate /home/users/eeslb/.bashrc
source activate /home/users/eeslb/miniconda3/envs/pyart_3_8_biodar

python plot_radar.py /gws/nopw/j04/ncas_radar_vol2/pestdar/india/raw_data/Delhi05/ 1 100 /gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/plots/india/delhi

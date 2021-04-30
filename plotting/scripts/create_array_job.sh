#!/bin/bash

job_directory=$PWD/
job_name=salalah_test
DATA_DIR=/gws/nopw/j04/ncas_radar_vol2/pestdar/oman/salalah/raw/cfradial/
OUT_DIR=/gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/plots/oman/salalah/

DATES=(${DATA_DIR}/*)

ARRAY_END=$((${#DATES[@]}-1))

JOB_FILE="${job_directory}/test_array.sbatch"

	echo "#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --job-name=${job_name}
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=5:00
#SBATCH --array=0-${ARRAY_END}

DATA_DIR=${DATA_DIR}

DATES=(\${DATA_DIR}/*)

DATE_PATH=\${DATES[\$SLURM_ARRAY_TASK_ID]}

source activate /home/users/eeslb/.bashrc
source activate /home/users/eeslb/miniconda3/envs/pyart_3_8_biodar

python plot_radar.py \${DATE_PATH} 1 100 ${OUT_DIR}" > $JOB_FILE

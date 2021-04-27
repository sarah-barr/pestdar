## Plotting

There are 3 ways to create plots depending on what you need:
1. Individual plots (single timestep) can be made using [Jupyter notebooks](https://github.com/sarah-barr/pestdar/tree/main/plotting/notebooks)
2. A single month (India) or day (Oman) can be plotted using individual jobs on JASMIN
3. Multiple months or days can be plotted using array jobs on JASMIN

### 1. Creating single plots using Jupyter notebooks
The ? notebook can be run in Jupyter to create individual plots and test changes to the scripts. The code is (almost) identical to the code in the .py executable files used to make plots on JASMIN so should be easy to transfer to those scripts. 

### 2. Plotting on JASMIN
The [plot_radar.py](https://github.com/sarah-barr/pestdar/blob/main/plotting/scripts/plot_radar.py) script is used to plot all files from a particular directory. The script takes 4 inputs:
```
plot_radar.py dpath ele r spath 
```
- dpath: path to directory containing raw radar data (ie a single month for India data or single day for Oman data)
- ele: elevation (index rather than actual value)
- r: radius around radar to plot (in km)
- spath: path to output directory. The script will create directories for the date/elevation/radius being plotted so this should be where you want those directories to be. 

Example:
```
plot_radar.py /gws/nopw/j04/ncas_radar_vol2/pestdar/india/raw_data/Delhi05/ 1 100 /gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/plots/india/delhi 
```

#### 2.1. Single jobs
The [run_single_job.sh](https://github.com/sarah-barr/pestdar/blob/main/scripts/run_single_job.sh) bash script is used to submit individual jobs to JASMIN and run the plot_radar.py script.
The job script and the python script should be in a directory together and then submitted from there using `sbatch run_single_job.sh`. The `.err` and `.out` output files will then be saved to the same directory. When the job is finished the `.out` file will tell you how many plots have been created:  
```
Plotting Delhi data: Delhi04
Directory created: /gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/plots/india/delhi/range_100km/elevation_2-00/Delhi04
4169 files found. 0 plots already existed. 4169 plots made. 0 errors. 4169 total plots
```
If these lines aren't there the job didn't finish for some reason and you should check in the `.err` file. The script will check if the plot already exists so if the job doesn't finish (eg it runs out of time) you can run it again without it making all the plots again. 

### 2.2. Job arrays
Job arrays are groups of jobs with the same executable and resource requirements, but different input files. They are used here to submit jobs simultaneously to plot multiple months/days. The 

The [create_array_job.sh](https://github.com/sarah-barr/pestdar/blob/main/scripts/create_array_job.sh) wrapper script can be used to create the job submission script. This takes the data directory (ie all salalah or all india data) and then checks how many directories there are within it and creates an array job submission script. `job_name`, `DATA_DIR`, `OUT_DIR` (lines 4,5,6) should be changed and then running it (`bash create_array_job.sh`) will write the [run_array_job.sbatch](https://github.com/sarah-barr/pestdar/blob/main/scripts/run_array_job.sbatch) script in the current working directory. This can then be submitted as above, the script will loop over the directories and submit each as a seperate task. Each task will have its own `.err` and `.out` output file. The job time limit relates to the time limit of each task. 

If you know how many jobs you need then you can just edit [run_array_job.sbatch](https://github.com/sarah-barr/pestdar/blob/main/scripts/run_array_job.sbatch) without using the wrapper. The array size (`#SBATCH --array=0-42`), the directory paths and possibly the time limit will need to be changed. 

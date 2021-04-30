# Radar gridding and writing geotiffs
Scripts here relate to mapping radar objects to a Cartesian grid and writing geotiffs.

### Gridding options 
There are multiple options for mapping radar objects to grids. The [test_grid_options](https://github.com/sarah-barr/pestdar/blob/main/geotiffs/notebooks/) notebook tests different options and settings. The best values from this testing were then used for creating the geotiffs

### Write geotiffs
Multiband geotiffs are created from grid objects using ?? script or ?? notebook. The .py executable script can be submitted to Jasmin using the run_geotiff_job.sh bash script. The script creates geotiffs for all data on a particular date, input and output paths are specified in the file. The output is a multiband geotiff where each band represents a different radar variable. This script requires gdal and rasterio to be installed. 


#!/usr/bin/env python
# coding: utf-8

'''
plot all radar data from one day 
inputs:
 date - date to plot 
 elevation - radar elevation to plot
 range - range from radar to plot (in km)

useage:
 plot_radar_all.py date elevation range

output:
 a plot is created for each timestep on the given date
 plots are saved ?? somewhere
 
'''
import os
import sys
import gc
import pyart as pyart
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from pyart.config import get_metadata
from glob import glob
import netCDF4
import cartopy
import cartopy.crs as ccrs
import numpy as np
import warnings
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from osgeo import gdal
from osgeo import osr
import rasterio
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
import matplotlib.colors as colors
import matplotlib.patches as mpatches

warnings.filterwarnings("ignore", category=DeprecationWarning) 


def get_args(args_in):
    """
    get input arguments
    """
    # check length of arguments:
    if len(args_in) != 5:
        # get name of the executable:
        self_name = os.path.basename(args_in[0])
        # print error and exit:
        sys.stderr.write('Error: incorrect or missing input variables')
        sys.exit()
    # get values:
    dpath = args_in[1]
    ele = int(args_in[2])
    r = int(args_in[3])
    spath = args_in[4]
    # return values:
    return dpath, ele, r, spath

def main():
    dpath, ele, r, spath = get_args(sys.argv)
    fnames = os.listdir(dpath)
    date = os.path.basename(os.path.normpath(dpath))
    
    nfiles = len(fnames)
    exists_count = 0
    plt_count = 0
    err_count = 0
    
    if 'Delhi' in dpath:
        print('Plotting Delhi data: {}'.format(date), flush = True)
        location = 'delhi'
        elevations = ([ 0.49987793,  0.99975586,  1.9995117 ,  2.9992676 ,  4.4989014 ,
                        5.998535  ,  8.997803  , 12.002563  , 16.001587  , 21.000366  ])
        elevation = elevations[ele]
        plt_params = delhi_plot_setup(os.path.join(dpath, fnames[0]), ele)
        spath = os.path.join(spath, location)
        
    elif 'salalah' in dpath:
        print('Plotting Salalah data: {}'.format(date), flush = True)
        elevations = [ 0.5,  1. ,  1.5,  2. ,  4. ,  6. , 10. , 15. ]
        elevation = elevations[ele]
        plt_params = salalah_plot_setup(os.path.join(dpath, fnames[0]), ele)
        

    elif 'JPR' in dpath:
        print('Plotting Jaipur data: {}'.format(date), flush = True)
        location = 'jaipur'
        elevations = ([ 0.49987793,  0.99975586,  1.9995117 ,  2.9992676 ,  4.4989014 ,
                        5.998535  ,  8.997803  , 12.002563  , 16.001587  , 21.000366  ])
        elevation = elevations[ele]
        plt_params = jpr_plot_setup(os.path.join(dpath, fnames[0]), ele)
        if not plt_params:
            err_count = err_count + 1
        else:
            spath = os.path.join(spath, location)
    else: 
        print('Unknown radar')
    
    for fname in fnames:
        fpath = os.path.join(dpath, fname)
        ele_str = '{:.2f}'.format(elevation)
        ele_str = ele_str.replace('.','-')
        sname = os.path.splitext(fname)[0] + '_{}km'.format(r) + '_ele{}'.format(ele_str) + '.png'
        sdir = os.path.join(spath, 'range_{}km'.format(r), 'elevation_{}'.format(ele_str), date)
        sfile = os.path.join(sdir, sname)
        if not os.path.exists(sdir):
            os.makedirs(sdir)
            print('Directory created: {}'.format(sdir))
            try:
                radar_plotting(fpath, ele, r, plt_params, sfile)
                plt_count += 1
            except:
                err_count += 1
                        
        else:
            if not os.path.isfile(sfile):
                try:
                    radar_plotting(fpath, ele, r, plt_params, sfile)
                    plt_count += 1
                except:
                    err_count += 1
            else:
                exists_count += 1
                continue
    tot_count = plt_count + exists_count
    print('{} files found. {} plots already existed. {} plots made. {} errors. {} total plots'. format(nfiles, exists_count, plt_count, err_count, tot_count), flush = True)

def jpr_plot_setup(fpath, ele):
    
    inc_var = ['differential_reflectivity',
             'reflectivity',
             'spectrum_width',
             'velocity']
    try:
        radar=pyart.io.read_sigmet(fpath, include_fields = inc_var)
        fields = sorted(list(radar.fields))

        interval = np.hstack([np.linspace(0, 0.45), np.linspace(0.55, 1)])
        cols = plt.cm.seismic(interval)
        cmap_vel = matplotlib.colors.LinearSegmentedColormap.from_list('name', cols)

        plt_names = ['Differential Reflectivity', 'Reflectivity', 'Spectrum Width', 'Doppler Velocity']

        plt_cols = ['pyart_RefDiff', 'pyart_HomeyerRainbow','pyart_NWSRef', cmap_vel]
        
        vmin = [-2, -30,0,-15]
        vmax = [10, 60,5,15]
        units = ['dB','dBZ','m/s','m/s']

        plt_order = [2,1,3,4]
        
        instrument_name = radar.metadata['instrument_name'].decode('utf-8')
        elevation = radar.fixed_angle['data'][ele]
        radar = []
        
        return inc_var, fields, plt_names, plt_cols, vmin, vmax, units, plt_order, instrument_name, elevation
    
    except OSError:
        err = 1 
        if err != 0:
            pass

def delhi_plot_setup(fpath, ele):
    '''
    load data and set plot names etc for delhi radar
    
    '''
    #specify variables to include (exclude_fields not working)
    inc_var = ['reflectivity', 
               'velocity', 
               'spectrum_width', 
               'differential_reflectivity', 
               'differential_phase', 
               'normalized_coherent_power', 
               'cross_correlation_ratio', 
               'specific_differential_phase']

    #read data
    radar=pyart.io.read(fpath, include_fields = inc_var)

    #define parameters for plotting
    fields = sorted(list(radar.fields))

    interval = np.hstack([np.linspace(0, 0.45), np.linspace(0.55, 1)])
    cols = plt.cm.seismic(interval)
    cmap_vel = matplotlib.colors.LinearSegmentedColormap.from_list('name', cols)

    plt_names = ['Cross Correlation Ratio',
                 'Differential Phase',
                 'Differential Reflectivity',
                 'Normalized Coherent Power (SQI)',
                 'Horizontal Reflectivity',
                 'Specific Differential Phase',
                 'Spectrum Width',
                 'Doppler Velocity',]

    plt_cols = ['pyart_RefDiff',
               'pyart_Wild25',
               'pyart_RefDiff',
               'pyart_Carbone17',
               'pyart_HomeyerRainbow',
               None,
              'pyart_NWSRef',
               cmap_vel]

    vmin = [0, 0,   -2, 0, -30, -1, 0, -15]
    vmax = [1, 180, 10, 1, 60, 7, 5, 15]
    units = ['', 'degrees', 'dB', '', 'dBZ', 'deg./km', 'm/s',  'm/s']
    plt_order = [3,7,2,4,1,8,5,6]

    instrument_name = radar.metadata['instrument_name'].decode('utf-8')
    elevation = radar.fixed_angle['data'][ele]
    radar = []
    
    return inc_var, fields, plt_names, plt_cols, vmin, vmax, units, plt_order, instrument_name, elevation
    
def salalah_plot_setup(fpath, ele):
    '''
    load data and set plot names etc for salalah radar
    
    '''
    #specify variables to include (exclude_fields not working)
    inc_var = ['RhoHV', 'SNR', 'SQI', 'V', 'W', 'ZDR', 'dBZ', 'uPhiDP']
    
    #read data
    radar=pyart.io.read(fpath, include_fields = inc_var)

    #define parameters for plotting
    fields = list(radar.fields)

    #change dop. vel. cmap 
    interval = np.hstack([np.linspace(0, 0.45), np.linspace(0.55, 1)])
    cols = plt.cm.seismic(interval)
    cmap_vel = matplotlib.colors.LinearSegmentedColormap.from_list('name', cols)
    
    plt_names = ['Cross Correlation Ratio',
                 'SNR',
                 'Normalized Coherent Power', 
                 'Doppler Velocity',
                 'Spectrum Width',
                 'Differential Reflectivity',
                 'Horizontal Reflectivity', 
                 'Unfiltered Differential Phase']
    
    plt_cols = ['pyart_RefDiff', 
                'pyart_Carbone17', 
                'pyart_Carbone17', 
                cmap_vel,
                'pyart_NWSRef', 
                'pyart_RefDiff',
                None, 
                'pyart_Wild25']
    
    units = ['','dB','','m/s','m/s','db','dBZ','degrees']
    vmin = [.5,0.,0.,-10,0,-2,-30,-360]
    vmax = [1.,100,1.,16.,5.,10.,60.,360]
    plt_order = [3,8,4,6,5,2,1,7]
    
    instrument_name = (radar.metadata['instrument_name'])[0:7]
    elevation = radar.fixed_angle['data'][ele]
    radar = []
    
    return inc_var, fields, plt_names, plt_cols, vmin, vmax, units, plt_order, instrument_name, elevation

def radar_plotting(fpath, ele, r, plt_params, sfile):
    
    inc_var, fields, plt_names, plt_cols, vmin, vmax, units, plt_order, instrument_name, elevation = plt_params
    
    radar=pyart.io.read(fpath, include_fields = inc_var)
    display = pyart.graph.RadarMapDisplay(radar)

    ## Plotting Options
    ele = ele #Elevation index
    r = r
    
    ## set the figure title
    instrument_name = instrument_name
    time_start = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    time_text = ' ' + time_start.strftime('%Y-%m-%d %H:%M:%SZ')
    elevation = elevation

    ## Figure Options
    
    nplts = len(fields)
    width=8 #in inches
    
    if nplts == 4:
        height = 2
        nrows = 1
    else:
        height=3.75 #in inches
        nrows = 2
        
    fig = plt.figure(figsize=(width, height), dpi = 200)
    ncols= int(nplts/nrows)

    ## Setting projection
    projection = ccrs.LambertConformal(central_latitude=radar.latitude['data'][0],
                                       central_longitude=radar.longitude['data'][0])

    ## adjust plot area (based on input range)
    #calculate lat lon from km 
    min_lat  = radar.latitude['data'][0]  - (r / 6378 ) * (180 / np.pi) #
    max_lat  = radar.latitude['data'][0]  + (r / 6378 ) * (180 / np.pi)
    min_lon = radar.longitude['data'][0] - (r / 6378 ) * (180 / np.pi) / np.cos(radar.latitude['data'][0] * np.pi/180)
    max_lon = radar.longitude['data'][0] + (r / 6378 ) * (180 / np.pi) / np.cos(radar.latitude['data'][0] * np.pi/180)

    ## get land cover info
    rst, img_ext, cmap, norm, name, lcp = get_land_cover(radar, projection)

    ## loop through fields to plot each one
    for field in fields:
        pnum = fields.index(field)

        ## create new axes for each field
        ax = plt.subplot(nrows,ncols,plt_order[pnum],projection=projection)
        plt.setp(ax.spines.values(), linewidth=0.25)

        ## add background image
        #ax.add_image(stamen_terrain, 8)

        ##format gridlines and labels
        xspace = 0.5
        yspace = 0.25
        lon_lines=np.arange(np.floor(min_lon), np.ceil(max_lon)+xspace, xspace)
        lat_lines=np.arange(np.floor(min_lat), np.ceil(max_lat)+yspace, yspace)
        gl = ax.gridlines(draw_labels=True, dms = True, x_inline=False, y_inline=False, alpha = 0.7, zorder = 0.5, linewidth = 0.25)
        gl.xlocator = mticker.FixedLocator(lon_lines)
        gl.ylocator = mticker.FixedLocator(lat_lines)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': '4'}
        gl.ylabel_style = {'size': '4'}
        gl.xpadding = 8

        ## some rubbish to turn the correct labels on/off
        if plt_order[pnum]==1:
            gl.right_labels = False
            #gl.bottom_labels = False
            gl.top_labels = False

        if plt_order[pnum]== 2 or plt_order[pnum]== 3:
            gl.left_labels = False
            gl.right_labels = False
            #gl.bottom_labels = False
            gl.top_labels = False

        if plt_order[pnum]== 6 or plt_order[pnum]== 7:
            gl.left_labels = False
            gl.right_labels = False
            gl.top_labels = False
            gl.bottom_labels = False

        if plt_order[pnum]== 4:
            gl.left_labels = False
            gl.right_labels = False
            #gl.bottom_labels = False
            gl.top_labels = False

        if plt_order[pnum]== 5:
            gl.right_labels = False
            gl.top_labels = False
            gl.bottom_labels = False

        if plt_order[pnum]== 8:
            gl.left_labels = False
            gl.right_labels = False
            gl.top_labels = False
            gl.bottom_labels = False

        ## plot each field
        display.plot_ppi_map(field,ele,ax=ax, vmin=vmin[pnum], vmax=vmax[pnum], title_flag= False, 
                             colorbar_flag = False, cmap = plt_cols[pnum] , embelish = False, alpha = 0.7,
                             min_lon = min_lon, max_lon = max_lon, min_lat = min_lat, max_lat = max_lat,
                             projection=projection,fig=fig, lat_0=radar.latitude['data'][0], lon_0=radar.longitude['data'][0])

        ## plot landcover raster
        ax.imshow(rst, extent = img_ext, cmap = cmap, norm =norm, interpolation='nearest', alpha = 0.5, zorder = 0.5)

        ## add title
        ax.text(0.03,0.97, plt_names[pnum], fontsize = 5, bbox=dict(facecolor=(1,1,1,0.9),
                                            edgecolor=(1,1,1,0.5),linewidth = 0.5, pad = 1.25), transform=ax.transAxes, va = 'top', ha = 'left', zorder = 10)


        ## format colorbar
        divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="7%", pad=0.075, axes_class = plt.Axes)

        if units[pnum]:
            height = '93%'
        else:
            height = '100%'

        cax =inset_axes(ax, width="7%",  # width = 5% of parent_bbox width
                       height=height,  # height : 50% 
                       loc='lower left',
                       bbox_to_anchor=(1.02, 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0)

        cbar = plt.colorbar(ax.collections[0], cax=cax, orientation="vertical")
        cbar.ax.tick_params(axis='y', direction='in', labelsize = 3.5, pad = 2, width = 0.25, length = 2)
        cbar.outline.set_linewidth(0.25)
        #cax.set_ylabel(radar.fields[field]['units'], rotation = 270, fontsize = 12,  weight="bold", labelpad = 4)
        ax.set_title(units[pnum], y = 0.94, x = 1.02, 
                      fontsize = 5, va = 'top', ha = 'left')

        ## add range rings
        range_rings = [25,50,75,100,125,150]
        plot_range_rings(range_rings, 45, radar, ax)

        ## add point for radar
        ax.plot(radar.longitude['data'][0], radar.latitude['data'][0], marker = 'o', markersize = 1, c = 'black', alpha = 0.7, transform = ccrs.PlateCarree())


    ## add main title 
    title = '$\\bf{Instrument\ name:}$' + '{}'.format(instrument_name) + '\n' + '$\\bf{Date\ and\ time:}$' + '{}'.format(time_text) + '\n' + '$\\bf{Elevation:}$' + '{:.2f}'.format(elevation)
    #title = 'Location: ' + '{}'.format(instrument_name) + '\n' + 'Data and Time: ' + '{}'.format(time_text) + '\n' + 'Elevation: ' + '{}'.format(ele)
    plt.suptitle(title, fontsize=5, y = 0.92, x = 0.13, ha = 'left', va = 'center')
    plt.subplots_adjust(hspace=0.001,wspace=0.25)   

    #land cover legend
    sort_inds = np.flip(np.argsort(lcp))
    patches = [mpatches.Patch(color=cmap(i), label='{} ({:.1f}%)'.format(name[i], lcp[i]), linewidth=0, alpha = 0.5) for i in sort_inds]
    fig.legend(handles=patches, fontsize = 4.5, loc='center right', bbox_to_anchor=(0.9, 0.92), borderaxespad=0., ncol=4, frameon=False)


    ## save figure
    sfile = sfile
    
    plt.savefig(sfile, format='png',bbox_inches='tight', dpi = 200)    
    plt.close('all')
    radar = []
    rst = []
    display = []
    gc.collect()
           
def plot_range_rings(range_rings_km, location, radar, ax):
    '''
    plot range rings at specified intervals 
    range_rings_km - list of rings to plot (in km)
    location - defines position of labels (0-360 degrees, 0 = top)
    '''
    loc_r = np.radians(location)
    transform = ccrs.AzimuthalEquidistant(central_longitude=radar.longitude['data'][0], 
                                          central_latitude=radar.latitude['data'][0])

    for ring in range_rings_km:
        angle = np.linspace(0., 2.0 * np.pi, 360)
        mask_angle = 22/ring #calculate angle needed for same length 'cut out' on each range ring (= arc length/r)

        for i in range(len(angle)):
            if loc_r-(mask_angle/2) <= angle[i] <= loc_r+(mask_angle/2):
                angle[i] =np.nan #nan values to exclude from plot where label will go

        label = (str(ring) + ' km')
        xpts = ring * 1000. * np.sin(angle)
        ypts = ring * 1000. * np.cos(angle)
        rot = 360 - location
        ax.plot(xpts, ypts, linestyle='--', c = 'gray', linewidth = 0.25, transform = transform, zorder = 10)
        txt = ax.text(ring*1000*np.sin(loc_r),ring*1000*np.cos(loc_r),label, va = 'center', ha = 'center', fontsize = 3, alpha = 0.7, 
                      color = 'gray', clip_on = True, transform = transform, rotation = rot)
        txt.clipbox = ax.bbox

def get_land_cover(radar, projection):
    
    lcpath = '/gws/nopw/j04/ncas_radar_vol1/eeslb/pestdar/data/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif'
    
    r = 100

    min_lat  = radar.latitude['data'][0]  - (r / 6378 ) * (180 / np.pi)
    max_lat  = radar.latitude['data'][0]  + (r / 6378 ) * (180 / np.pi)
    min_lon = radar.longitude['data'][0] - (r / 6378 ) * (180 / np.pi) / np.cos(radar.latitude['data'][0] * np.pi/180)
    max_lon = radar.longitude['data'][0] + (r / 6378 ) * (180 / np.pi) / np.cos(radar.latitude['data'][0] * np.pi/180)

    #get data from window based on map extent 
    
    with rasterio.open(lcpath) as src:
        src_cmap = src.colormap(1)
        rst = src.read(1, window=from_bounds(min_lon, min_lat, max_lon, max_lat, src.transform))
        gt = src.window_transform(from_bounds(min_lon, min_lat, max_lon, max_lat, src.transform))

    #reclassify forest subtypes
    mask = np.where(np.logical_and(rst>=111, rst<=126))
    rst[mask] = 10

    unique_vals = np.unique(rst)
    
    #group low % types
    for i in unique_vals:
        total_count = np.size(rst)
        val_count = np.size(rst[rst == i])
        per = (val_count/total_count)*100
        if per < 0.1:
            rst[rst == i] = 0

    unique_vals = np.unique(rst)

    class_names = ['Other', 'Forest', 'Shrubland', 'Herbaceous veg.', 'Cropland', 'Built-up', 'Bare/sparse veg.', 'Snow & ice', 
                   'Water', 'Wetland', 'Moss/lichen', 'Open sea']

    class_values = ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100', '200']
    col = ['lightgray', 'darkgreen', 'brown', 'olive', 'wheat', 'black', 'blanchedalmond', 'white', 'lightblue', 'cadetblue', 'white', 'lightblue']


    #get names and colors
    name = []
    c = []
    lcp = []
    for i in range(len(unique_vals)):
        data_ind = unique_vals[i]
        ind = class_values.index(str(data_ind))
        name.append(class_names[ind])
        c.append(col[ind])
        new_val = len(c)-1
        rst[rst == data_ind] = new_val
        lcp.append((np.size(rst[rst == new_val])/np.size(rst))*100)
        #print('Name: ', class_names[ind], ' Colour: ', col[ind], new_val)

    #create colormap and set bounds 
    cmap = matplotlib.colors.ListedColormap(c)
    bounds = np.unique(rst)
    norm = colors.BoundaryNorm(np.arange(len(bounds)+1)-0.5,cmap.N)

    ## reproject bounds to lambert conformal and set raster extent
    lower_left = projection.transform_point(min_lon, min_lat, ccrs.PlateCarree())
    upper_right = projection.transform_point(max_lon, max_lat, ccrs.PlateCarree())

    img_ext = (lower_left[0], upper_right[0], lower_left[1], upper_right[0])

    ## plot
   # plot = ax.imshow(rst, extent = img_ext, cmap = cmap, norm = norm, interpolation='nearest', alpha = 0.7)
   # ax.gridlines(ccrs.PlateCarree(), draw_labels = True)

    ## add legend 
    #patches = [mpatches.Patch(color=cmap(i), label=name[i], alpha = 0.7) for i in range(len(name))]
    # put those patched as legend-handles into the legend
    #plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )

    return rst, img_ext, cmap, norm, name, lcp
            
if __name__ == '__main__':
    main()

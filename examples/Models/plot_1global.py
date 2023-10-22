#!/usr/bin/env python
"""This module contains an example of plotting a map using avni codes"""

import os
import argparse #parsing arguments
import matplotlib.pyplot as plt
import ntpath
import xarray as xr
import numpy as np

########################### IMPORT AVNI MODULES   #####################################

from avni.tools import stage,get_fullpath,get_filedir
from avni.data import update_file
from avni.models import readepixfile
from avni.plots import globalmap

#########################################################
def main():
    parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
    parser.add_argument('-m', '--file', type=str, default=
        'S362ANI+M.BOX25km_PIX1X1.avni.nc4',help='Model file in nc4 or epix')
    parser.add_argument('-p', '--parameter', type=str, default='vs',
        help='Parameter of interest')
    parser.add_argument('-d', '--depth', type=float, default=150.0,
        help='Depth request in km')
    parser.add_argument('-j', '--projection', type=str,default='robin',
        help='map projection')
    parser.add_argument('-u', '--upper_bound', type=float, default=6.0,
        help='Upper bound for color scale saturation level (percent)')
    parser.add_argument('-l', '--lower_bound', type=float, default=-6.0,
        help='Lower bound for color scale saturation level (percent)')
    parser.add_argument('-c', '--color', type=str, default='avni',
        help='Color palette e.g. avni or bk')
    parser.add_argument('-r', '--resolution', type=str, default='l',
        help='Resolution')
    parser.add_argument('-t', '--colorcontour', type=int, default=20,
        help='color contour levels')
    parser.add_argument('-f', '--format', type=str,default='png',
        help='Output file format')
    parser.add_argument('-o', '--output', action='store_true',
        help='Save the figures in files')
    parser.add_argument('-s', '--shading', action='store_true',
        help='Shade topography')
    arg = parser.parse_args()

    try:
        # stage the file for plotting
        if not os.path.isfile(get_fullpath(arg.file)):
            raise IOError('File does not exist locally :'+arg.file)
        ierror = stage(get_fullpath(arg.file),overwrite=True)
    except:
        # update the file from the server
        localfile, success = update_file(arg.file,subdirectory='MODELS/S362ANI+M')
        if success: print('....Downloaded file from AVNI server as '+localfile)

    # Read the file
    model = ntpath.basename(arg.file)
    if not os.path.isfile(get_filedir()+'/'+model):
        raise IOError('File does not exist locally: '+get_filedir()+'/'+model)
    try:
        latlonval,metadata,_ = readepixfile(get_filedir()+'/'+model)
        try:
            title = 'Depth : '+metadata['DEPTH_RANGE']+' km'
            outfile = model+'.'+arg.parameter+'.'+metadata['DEPTH_IN_KM']+'km.'+arg.format
        except:
            title = arg.parameter
            outfile = model+'.'+arg.parameter+'.'+arg.format
        metadata['WHAT'] = arg.parameter
        metadata['UNIT'] = ''
    except:
        metadata = {}
        ds = xr.open_dataset(get_filedir()+'/'+model)
        metadata['WHAT'] = arg.parameter
        metadata['UNIT'] = '%' if ds[arg.parameter].attrs['unit']=='percent' else ds['vs'].attrs['unit']
        outfile = model.split('avni.nc4')[0]+arg.parameter+'.'+str(arg.depth)+'km.'+arg.format

        # find the index, choosing the first if multiple ones available
        start_depths = ds[arg.parameter].attrs['start_depths']
        end_depths = ds[arg.parameter].attrs['end_depths']
        find = np.logical_and(start_depths <= arg.depth,end_depths >= arg.depth)
        if not np.any(find): raise ValueError('depth queried '+str(arg.depth)+' not found among the limits : ',np.vstack((start_depths,end_depths)).T)
        depindex = np.where(find)[0][0]
        start_depth = start_depths[depindex]
        end_depth = end_depths[depindex]
        title = 'Depth : '+str(arg.depth)+' km ['+str(start_depth)+' - '+str(end_depth)+' km]'
        latlonval = ds[arg.parameter][depindex]

    # Plot the file
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    projection=arg.projection; vmin = arg.lower_bound; vmax = arg.upper_bound
    if projection=='ortho':
        globalmap(ax,latlonval,vmin,vmax,grid=[30.,30.],gridwidth=1,projection=projection,colorlabel=metadata['WHAT']+' ('+metadata['UNIT']+ ')',lat_0=0,lon_0=150,colorcontour=arg.colorcontour,colorpalette=arg.color,shading=arg.shading,resolution=arg.resolution)
    elif projection=='robin':
        globalmap(ax,latlonval,vmin,vmax,grid=[30.,90.],gridwidth=0,projection=projection,colorlabel=metadata['WHAT']+' ('+metadata['UNIT']+ ')',lon_0=150,colorcontour=arg.colorcontour,colorpalette=arg.color,shading=arg.shading,resolution=arg.resolution)
    else:
        globalmap(ax,latlonval,vmin,vmax,grid=[30.,90.],gridwidth=0,projection=projection,colorlabel=metadata['WHAT']+' ('+metadata['UNIT']+ ')',lat_0=0,lon_0=150,colorcontour=arg.colorcontour,colorpalette=arg.color,shading=arg.shading,resolution=arg.resolution)
    ax.set_title(title)
    if arg.output:
        fig.savefig(outfile,dpi=300)
    else:
        plt.show()
    return

if __name__== "__main__":
    main()

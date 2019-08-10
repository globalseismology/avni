#!/usr/bin/env python
"""This module contains an example of plotting a map using rem3d codes"""

import argparse #parsing arguments
import matplotlib.pyplot as plt
import ntpath
import xarray as xr
import pdb

########################### IMPORT REM3D MODULES   #####################################
import rem3d
#########################################################
def main():
    parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
    parser.add_argument('-m', '--file', type=str, default=
        'S362ANI+M.BOX25km_PIX1X1.rem3d.nc4',help='Model file in nc4 or epix')
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
    parser.add_argument('-c', '--color', type=str, default='rem3d',
        help='Color palette e.g. rem3d or bk')
    parser.add_argument('-t', '--colorcontour', type=int, default=20,
        help='color contour levels')
    parser.add_argument('-o', '--format', type=str,default='png',
        help='Output file format')
    arg = parser.parse_args()

    try:
        # stage the file for plotting
        rem3d.tools.stage(arg.file)
    except:
        # update the file from the server
        rem3d.data.update_file(arg.file)
    model = ntpath.basename(arg.file)

    # Read the file
    try:
        latlonval,metadata,_ = rem3d.models.readepixfile(rem3d.tools.get_filedir()+'/'+model)
        title = 'Depth : '+metadata['DEPTH_RANGE']+' km'
        outfile = model+'.'+arg.parameter+'.'+metadata['DEPTH_IN_KM']+'km.'+arg.format
    except:
        raise NotImplementedError('netcdf not implemented yet')
        metadata = {}
        ds = xr.open_dataset(rem3d.tools.get_filedir()+'/'+model)
        metadata['WHAT'] = arg.parameter
        metadata['UNIT'] = '%' if ds[arg.parameter].attrs['unit']=='percent' else ds['vs'].attrs['unit']
        outfile = model+'.'+arg.parameter+'.'+str(arg.depth)+'km.'+arg.format
        title = 'Depth : '+str(arg.depth)+' km'
        #title = 'Depth : '+str(arg.depth)+' km ['+str+' - '+str+' km]'

    # Plot the file
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    projection=arg.projection; vmin = arg.lower_bound; vmax = arg.upper_bound
    if projection=='ortho':
        rem3d.plots.globalmap(ax,latlonval,vmin,vmax,grid=[30.,30.],gridwidth=1,projection=projection,colorlabel=metadata['WHAT']+' ('+metadata['UNIT']+ ')',lat_0=0,lon_0=150,colorcontour=arg.colorcontour,colorpalette=arg.color)
    elif projection=='robin':
        rem3d.plots.globalmap(ax,latlonval,vmin,vmax,grid=[30.,90.],gridwidth=0,projection=projection,colorlabel=metadata['WHAT']+' ('+metadata['UNIT']+ ')',lon_0=150,colorcontour=arg.colorcontour,colorpalette=arg.color)
    else:
        rem3d.plots.globalmap(ax,latlonval,vmin,vmax,grid=[30.,90.],gridwidth=0,projection=projection,colorlabel=metadata['WHAT']+' ('+metadata['UNIT']+ ')',lat_0=0,lon_0=150,colorcontour=arg.colorcontour,colorpalette=arg.color)
    ax.set_title(title)
    plt.show()
    fig.savefig(outfile,dpi=300)

    return

if __name__== "__main__":
    main()

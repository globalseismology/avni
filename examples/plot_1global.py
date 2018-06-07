#!/usr/bin/env python
"""This module contains an example of plotting a map using rem3d codes"""

import sys,os
import argparse #parsing arguments
import matplotlib.pyplot as plt

################################ IMPORT REM3D MODULES   #####################################
import rem3d

#########################################################
def main():
    parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
    parser.add_argument('-m', '--model', type=str,default='ME16',
        help='3D model')
    parser.add_argument('-j', '--projection', type=str,default='robin',
        help='map projection')
    parser.add_argument('-u', '--upper_bound', type=float, default=4.0,
        help='Upper bound for color scale saturation level (percent)')
    parser.add_argument('-l', '--lower_bound', type=float, default=-4.0,
        help='Lower bound for color scale saturation level (percent)')
    parser.add_argument('-p', '--parameter', type=str, default='Vp',
        help='Parameter to plot - Vs, Vp, etc.')
    parser.add_argument('-d', '--depth', type=float, default=100,
        help='Depth in km')
    parser.add_argument('-o', '--format', type=str,default='.pdf',
        help='Output file format')
    parser.add_argument('-b', '--base_path', type=str, default='.',
        help='Database path containing all files')
    arg = parser.parse_args()
    
    filename=arg.model+'.'+arg.parameter+'.'+str(arg.depth)+'km.epix'
         
    # Read the file
    ierror, latlonval=rem3d.models.readepixfile(arg.base_path+'/'+filename)
    
    # Plot the file
    fig=plt.figure() 
    ax=fig.add_subplot(1,1,1)
    projection=arg.projection; vmin = arg.lower_bound; vmax = arg.upper_bound
    if projection=='ortho':
        rem3d.plots.globalmap(ax,latlonval,vmin,vmax,grid=[30.,30.],gridwidth=1,projection=projection,colorlabel=arg.parameter+' (%)',lat_0=0,lon_0=150,colorpalette='rem3d')
    else:
        rem3d.plots.globalmap(ax,latlonval,vmin,vmax,grid=[30.,90.],gridwidth=0,projection=projection,colorlabel=arg.parameter+' (%)',lat_0=0,lon_0=150,colorpalette='rem3d')
    ax.set_title('Depth : '+str(arg.depth)+' km')
    plt.show()
    fig.savefig(filename+arg.format,dpi=300)
        
    return

if __name__== "__main__":
    main()

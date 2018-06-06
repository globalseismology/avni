#!/usr/bin/env python
"""This module contains the various subroutines used for plotting
Usage import """
#####################  IMPORT STANDARD MODULES   ######################################   

from __future__ import division
from math import cos, pi, log, sin, tan, atan, atan2, sqrt, radians, degrees, asin, modf
import sys,os
import argparse #parsing arguments
import numpy as np #for numerical analysis
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import multiprocessing
import cartopy.crs as ccrs
from joblib import Parallel, delayed
import pdb    #for the debugger pdb.set_trace()
# from scipy.io import netcdf_file as netcdf #reading netcdf files
from netCDF4 import Dataset as netcdf #reading netcdf files
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import time
import progressbar
# For polar sectionplot
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator,DictFormatter,FixedLocator
from matplotlib import gridspec # Relative size of subplots

####################       IMPORT OWN MODULES     ######################################
from . import mapping
from . import tools
from . import data
from . import models
from . import constants
########################      GENERIC   ################################################                       
                    
def get_colors(val,xmin=-1.,xmax=1.,palette='coolwarm'):
    """gets the value of color for a given palette"""
    jet = cm = get_cmap(palette) 
    cNorm  = mcolors.Normalize(vmin=xmin, vmax=xmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=palette)
    colorVal = scalarMap.to_rgba(val)
    return colorVal
    
def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = cm = get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))
    
    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]
    
    return cmap.from_list(cmap.name + "_gray", colors, cmap.N)

def make_colormap(seq,name='CustomMap'):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap(name, cdict)

def getcolorlist(cptfile):
    """Get a tuple for colorlist from a cptfile"""
    currentdir=os.getcwd()
    try: 
        f = open(cptfile, 'r')
    except IOError:
        print "File ("+cptfile+") does not exist in the current directory - "+currentdir
        sys.exit(2)
    cptarr=np.genfromtxt(cptfile, dtype=None,comments="#")
    colorlist=[]
    for irow in np.arange(len(cptarr)): 
        tups=cptarr[irow][1]/255.,cptarr[irow][2]/255.,cptarr[irow][3]/255.
        val=(cptarr[irow][4]-cptarr[0][4])/(cptarr[len(cptarr)-1][0]-cptarr[0][4])
        if irow==1:
            colorlist.append(tups)
        elif irow > 1 and irow < len(cptarr)-1:
            colorlist.append(tups)    
            colorlist.append(val)    
            colorlist.append(tups)    
    return colorlist
    
def customcolorpalette(name='rem3d',cptfolder='~/CPT',colorlist=None,colormax=2.,middlelimit=0.5,ifgraytest=0):
    """Used to return preset color palettes from cptfolder. ifgraytest test how the figure looks in gray scale. (-colormax,colormax) are the limits of the colorbar. zerolimit is the limit to which the middle color (e.g. grey) will extend on either side of colorttmax mid. """
    c = mcolors.ColorConverter().to_rgb    
    if name=='r_lgrey_b':
        colorlist=[c('blue'), c('lightgray'), (2.*colormax-2.*middlelimit)/(4.*colormax), c('lightgray'),c('lightgray'), (2.*colormax+2.*middlelimit)/(4.*colormax), c('lightgray'),c('red'), 1., c('red')]
    elif name=='rem3d':
        cptfolder=tools.get_fullpath(cptfolder)
        if os.path.isfile(cptfolder+'/bk1_0.cpt_'):
            colorlist=getcolorlist(cptfolder+'/bk1_0.cpt_')
        else:
            print "Error: Could not find file "+cptfolder+'/bk1_0.cpt_'    
            pdb.set_trace()
            sys.exit(2)
    elif name=='hit1':
        cptfolder=tools.get_fullpath(cptfolder)
        if os.path.isfile(cptfolder+'/hit1.cpt_'):
            colorlist=getcolorlist(cptfolder+'/hit1.cpt_')
        else:
            print "Error: Could not find file "+cptfolder+'/hit1.cpt_'    
            pdb.set_trace()
            sys.exit(2)
    elif name=='yuguinv':
        cptfolder=tools.get_fullpath(cptfolder)
        if os.path.isfile(cptfolder+'/yu1_2inv.new.cpt_'):
            colorlist=getcolorlist(cptfolder+'/yu1_2inv.new.cpt_')
        else:
            print "Error: Could not find file "+cptfolder+'/yu1_2inv.new.cpt_'
            pdb.set_trace()
            sys.exit(2)
        
    if colorlist is None: sys.exit("No colorlist found")
    custom_cmap = make_colormap(colorlist,name)
    cmx.register_cmap(name=custom_cmap.name, cmap=custom_cmap)
    palette=custom_cmap.name
    
    if ifgraytest==1:
        palette=grayify_cmap(palette)
        
    return custom_cmap    

def standardcolorpalette(name='rem3d',RGBoption='rem3d',reverse=True):
    """
    Get a custom REM3D color palette from constants.py
    """
    if reverse:
        RGBlist=constants.colorscale[RGBoption]['RGB'][::-1]
    else:
        RGBlist=constants.colorscale[RGBoption]['RGB']
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(name, RGBlist,N=len(RGBlist))
    cmx.register_cmap(name=custom_cmap.name, cmap=custom_cmap)
    return custom_cmap    

                
############################### PLOTTING ROUTINES ################################        

def plot_hotspots(m, dbs_path = '.', lon360 = False, **kwargs):
    """Reads hotspots.pkl from dbs_path and plots on to map index m lon360 is False if the no lon above 180 is permitted. **kwargs denotes the arguments, if any, for scatter."""

    # Earlier was in pickle format, cross-platform compatibility required json
    # hotspots = pickle.load(open('%s/hotspots.pkl' % (dbs_path), 'rb'))
    # tools.writejson(hotspots,'%s/hotspots.json' % (dbs_path))
    hotspots = tools.readjson('%s/hotspots.json' % (dbs_path))
    if lon360:
        hotspots[:,0] = (hotspots[:,0] + 360) % 360.0
    x, y = m(hotspots[:,0], hotspots[:,1])
    if kwargs:
        m.scatter(x, y, **kwargs)
    else:
        m.scatter(x, y)
    return    

def plot_plates(m, dbs_path = '.', lon360 = False, boundtypes=['ridge', 'transform', 'trench'],**kwargs):
    """Reads hotspots.pkl from base_path and plots on to map index m lon360 is False if the no lon above 180 is permitted. **kwargs denotes the arguments, if any, for scatter."""
    
    # Earlier was in pickle format, cross-platform compatibility required json
    # ridge,ridgeloc=pickle.load(open('%s/ridge.pkl' % (dbs_path),'rb'))
    # tools.writejson(np.array([ridge,ridgeloc.tolist()]),'%s/ridge.json' % (dbs_path))
    for bound in boundtypes:
        #name, segs = pickle.load(open('%s/%s.pkl' % (dbs_path,bound), 'rb'))
        name, segs = tools.readjson('%s/%s.json' % (dbs_path,bound))
        segs=np.array(segs)
        ind_nan, = np.nonzero(np.isnan(segs[:,0]))
        segs[ind_nan,0] = 0
        segs[ind_nan,1] = 0
        if lon360:
            segs[:,0] = (segs[:,0] + 360) % 360.0
        x, y = m(segs[:,0], segs[:,1])
        x[ind_nan] = np.nan
        y[ind_nan] = np.nan
        dx = np.abs(x[1:] - x[:-1])
        ind_jump, = np.nonzero(dx > 1000000)
        x[ind_jump] = np.nan
        if kwargs:
            m.plot(x, y, '-', **kwargs)
        else:
            m.plot(x, y, '-')
    return

def globalmap(ax,valarray,vmin,vmax,dbs_path='.',colorlabel=None,colorpalette='rem3d',colorcontour=21,hotspots=False,grid=[30.,90.],gridwidth=0, **kwargs):
    """
    Plots a 2-D cross-section of a 3D model on a predefined axis ax.
    
    Parameters
    ----------
    
    latlonval : a named numpy array containing latitudes (lat), longitudes (lon) 
                and values (val). Can be initialized from three numpy arrays lat, lon and val
                $ data = np.vstack((lat,lon,val)).transpose()
                $ dt = {'names':['lat', 'lon', 'val'], 'formats':[np.float, np.float, np.float]}
                $ latlonval = np.zeros(len(data), dtype=dt)
                $ latlonval['lat'] = data[:,0]; latlonval['lon'] = data[:,1]; latlonval['val'] = data[:,2]
    
    vmin, vmax : minimum and maximum value of the color scale
    
    dbs_path : database path containing hotpot locations, coastlines etc.
    
    colorpalette : matploblib color scales or the REM3D one (default)
    
    colorcontour :  the number of contours for colors in the plot. Maximum is 520 and odd values
                    are preferred so that mid value is at white/yellow or other neutral colors.
    
    projection : map projection for the global plot
    
    colorlabel : label to use for the colorbar
    
    lat_0, lon_0 : center latitude and longitude for the plot
    
    outformat : format of the output file 
    
    kwargs : optional arguments for Basemap 
    """ 

    # set up map
    if kwargs:
        m = Basemap(ax=ax, **kwargs)
    else:
        m = Basemap(ax=ax,projection='robin', lat_0=0, lon_0=150, resolution='c')
    clip_path = m.drawmapboundary()
    m.drawcoastlines(linewidth=1.5)
    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(-90.,90.,grid[0])
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[True,False,False,False],linewidth=gridwidth,fontsize=14)
    meridians = np.arange(-180.,180.,grid[1])
    m.drawmeridians(meridians,labels=[False,False,False,True],linewidth=gridwidth,fontsize=14)

    # plot the model
    for ii in np.arange(len(valarray['lon'])): 
        if valarray['lon'][ii] > 180.: valarray['lon'][ii]=valarray['lon'][ii]-360.
    numlon=len(np.unique(valarray['lon']))
    numlat=len(np.unique(valarray['lat']))
    # grid spacing assuming a even grid
    grid_spacing = min(max(np.ediff1d(valarray['lat'])),max(np.ediff1d(valarray['lon'])))
    lat = np.arange(-90.+grid_spacing/2.,90.+grid_spacing/2.,grid_spacing)
    lon = np.arange(-180.+grid_spacing/2.,180.+grid_spacing/2.,grid_spacing)
    X,Y=np.meshgrid(lon,lat)
    val = np.empty_like(X)
    val[:] = np.nan;
    for i in xrange(0, valarray['lat'].size):
        ilon = np.where(X[0,:]==valarray['lon'][i])[0][0]
        ilat = np.where(Y[:,0]==valarray['lat'][i])[0][0]
        val[ilat,ilon] = valarray['val'][i]
    s = m.transform_scalar(val,lon,lat, 1000, 500)
    # Get the color map
    try:
        cpalette = plt.get_cmap(colorpalette)
    except ValueError:
        try:
            cpalette=standardcolorpalette(RGBoption=colorpalette)
        except KeyError:
            cpalette=customcolorpalette(colorpalette)
    
    # define the 10 bins and normalize
    if isinstance(colorcontour,np.ndarray) or isinstance(colorcontour,list): # A list of boundaries for color bar
        if isinstance(colorcontour,list): 
            bounds = np.array(colorcontour)
        else:
            bounds = colorcontour
        bounds = bounds[np.where(bounds < vmax)]
        mytks = np.append(bounds[bounds.nonzero()],np.ceil(vmax))
        bounds = np.append(bounds,np.ceil(vmax))
        spacing='uniform'
    elif isinstance(colorcontour,int): # Number of intervals for color bar
        bounds = np.linspace(vmin,vmax,colorcontour+1)
        mytks = np.arange(vmin,vmax+(vmax-vmin)/4.,(vmax-vmin)/4.)
        mytkslabel = [str(a) for a in mytks]
        spacing='proportional'
    else:
        print "Error: Undefined colorcontour, should be a numpy array, list or integer "
        sys.exit(2)        
    norm = mcolors.BoundaryNorm(bounds,cpalette.N)
    im = m.imshow(s, cmap=cpalette.name, clip_path=clip_path, vmin=vmin, vmax=vmax, norm=norm)
#    # add plates and hotspots
    dbs_path=tools.get_fullpath(dbs_path)
    plot_plates(m, dbs_path=dbs_path, color='w', linewidth=1.5)
    if hotspots: plot_hotspots(m, dbs_path=dbs_path, s=30, color='m', edgecolor='k')

#   Add a colorbar
    if colorlabel is not None:
#         cb = plt.colorbar(im,orientation='vertical',fraction=0.05,pad=0.05)
#         cb.set_label(colorlabel)
        # Set colorbar, aspect ratio
        cbar = plt.colorbar(im, alpha=0.05, aspect=12, shrink=0.5,norm=norm, spacing=spacing, ticks=bounds, boundaries=bounds,extendrect=True)
        cbar.solids.set_edgecolor("face")
        # Remove colorbar container frame
#         cbar.outline.set_visible(False)
        # Fontsize for colorbar ticklabels
        cbar.ax.tick_params(labelsize=14)
        # Customize colorbar tick labels
        cbar.set_ticks(mytks)
        mytkslabels = [str(int(a)) if (a).is_integer() else str(a) for a in mytks]
        cbar.ax.set_yticklabels(mytkslabels)
        # Colorbar label, customize fontsize and distance to colorbar
        cbar.set_label(colorlabel,rotation=90, fontsize=16, labelpad=5)
        # Remove color bar tick lines, while keeping the tick labels
#         cbarytks = plt.getp(cbar.ax.axes, 'yticklines')
#         plt.setp(cbarytks, visible=False)
        
    m.drawmapboundary(linewidth=1.5)    
    return m    

def plot1globalmap(latlonval,vmin,vmax,dbs_path='.',colorpalette='rem3d',projection='robin',colorlabel="Anomaly (%)",lat_0=0,lon_0=150,outformat='.pdf',ifshow=False):
    """
    Plot one global map and write it to a file or display on screen.
    
    Parameters
    ----------
    
    latlonval : a named numpy array containing latitudes (lat), longitudes (lon) 
                and values (val). Can be initialized from three numpy arrays lat, lon and val
                $ data = np.vstack((lat,lon,val)).transpose()
                $ dt = {'names':['lat', 'lon', 'val'], 'formats':[np.float, np.float, np.float]}
                $ latlonval = np.zeros(len(data), dtype=dt)
                $ latlonval['lat'] = data[:,0]; latlonval['lon'] = data[:,1]; latlonval['val'] = data[:,2]
    
    vmin, vmax : minimum and maximum value of the color scale
    
    dbs_path : database path containing hotpot locations, coastlines etc.
    
    colorpalette : matploblib color scales or the REM3D one (default)
    
    projection : map projection for the global plot
    
    colorlabel : label to use for the colorbar
    
    lat_0, lon_0 : center latitude and longitude for the plot
    
    outformat : format of the output file 
    
    ifshow : display the plot to the user if True
    """ 
    fig=plt.figure() 
    ax=fig.add_subplot(1,1,1)
    epixarr=models.readepixfile(filename)
    if projection=='ortho':
        globalmap(ax,epixarr,vmin,vmax,dbs_path,colorlabel,grid=[30.,30.],gridwidth=1,projection=projection,lat_0=lat_0, lon_0=lon_0,colorpalette=colorpalette)
    else:
        globalmap(ax,epixarr,vmin,vmax,dbs_path,colorlabel,grid=[30.,90.],gridwidth=0,projection=projection,lat_0=lat_0, lon_0=lon_0,colorpalette=colorpalette)
    if ifshow: plt.show()
    fig.savefig(filename+outformat,dpi=300)
    return 

    
    
    

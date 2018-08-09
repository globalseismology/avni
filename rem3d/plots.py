#!/usr/bin/env python
"""
This module contains the various subroutines used for plotting
Usage import
"""
#####################  IMPORT STANDARD MODULES   ######################################   

from __future__ import division
from math import cos, pi, log, sin, tan, atan, atan2, sqrt, radians, degrees, asin, modf
import sys,os
import numpy as np #for numerical analysis
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import multiprocessing
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
from . import tools
from . import data
from . import models
from . import constants
from . import mapping
########################      GENERIC   ################################################                       

def updatefont(ax=None,fontsize=15,fontname='sans-serif'): 
    """
    Updates the font type and sizes globally or for a particular axis handle
    
    Parameters
    ----------
    
    ax :  figure axis handle
    
    fontsize,fontname : font parameters
    
    Return:
    ----------
    
    ax : updated axis handle if ax is not None
    
    """
    if ax is None:
        plt.rcParams["font.family"] = fontname
        plt.rcParams["font.size"] = fontsize
    else:
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
            item.set_fontname(fontname)
        for item in ([ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(fontsize+2)
            item.set_fontname(fontname)
        ax.title.set_fontsize(fontsize+3)
        ax.title.set_fontname(fontname)
    return ax if ax is not None else None
                    
def standardcolorpalette(name='rem3d',RGBoption='rem3d',reverse=True):
    """
    Get a custom REM3D color palette from constants.py
    
    Parameters
    ----------
    
    name : color palette name that will be used elsewhere
    
    RGBoption : option for values of RGB from constants
    
    reverse: if the colors need to be reversed from those provided.
    
    """
    if reverse:
        RGBlist=constants.colorscale[RGBoption]['RGB'][::-1]
    else:
        RGBlist=constants.colorscale[RGBoption]['RGB']
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(name, RGBlist,N=len(RGBlist))
    cmx.register_cmap(name=custom_cmap.name, cmap=custom_cmap)
    return custom_cmap    
    
def get_colors(val,xmin=-1.,xmax=1.,palette='coolwarm',colorcontour=20):
    """gets the value of color for a given palette"""
    jet = cm = cmx.get_cmap(palette) 
    #cNorm  = mcolors.Normalize(vmin=xmin, vmax=xmax)
    bounds = np.linspace(xmin,xmax,colorcontour+1)
    cNorm = mcolors.BoundaryNorm(bounds,cm.N)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=palette)
    colorVal = scalarMap.to_rgba(val)
    return colorVal
                
############################### PLOTTING ROUTINES ################################        

def plot_hotspots(m, dbs_path = '.', lon360 = False, **kwargs):
    """
    Reads hotspots.pkl from dbs_path and plots on to map index m 

    Earlier, the data was in pickle format, cross-platform compatibility required json
    # hotspots = pickle.load(open('%s/hotspots.pkl' % (dbs_path), 'rb'))
    # tools.writejson(hotspots,'%s/hotspots.json' % (dbs_path))
    
    Parameters
    ----------
    
    m : figure axis handle
    
    dbs_path: path specified by user where hotspots.json is located. If not found,
              defaults to downloading the file from the  REM3D server.
              
    lon360 : is False if the no longitude above 180 is permitted and is wrapped around. 
    
    kwargs : denotes the arguments, if any, for scatter
               
    """
    
    try:
        hotspots = tools.readjson('%s/hotspots.json' % (dbs_path))
    except IOError: #Download to default directory
        filedir = tools.get_filedir(checkwrite=True,makedir=True)
        data.update_file('hotspots.json')
        hotspots = tools.readjson('%s/hotspots.json' % (filedir))
       
    if lon360:
        hotspots[:,0] = (hotspots[:,0] + 360) % 360.0
    x, y = m(hotspots[:,0], hotspots[:,1])
    if kwargs:
        m.scatter(x, y, **kwargs)
    else:
        m.scatter(x, y)
    return    

def plot_plates(m, dbs_path = '.', lon360 = False, boundtypes=['ridge', 'transform', 'trench'],**kwargs):
    """
    Plots different types of tectonic plates
    
    Parameters
    ----------
    
    m : figure axis handle
    
    dbs_path : path specified by user where hotspots.json is located. If not found,
              defaults to downloading the file from the  REM3D server.
              
    boundtypes : plate boundary types that will be plotted. Default are ridge, transform
                 and trench 
    
    lon360 : is False if the no longitude above 180 is permitted and is wrapped around. 
    
    kwargs : denotes the arguments, if any, for scatter
    
    """
    
    # Earlier was in pickle format, cross-platform compatibility required json
    # ridge,ridgeloc=pickle.load(open('%s/ridge.pkl' % (dbs_path),'rb'))
    # tools.writejson(np.array([ridge,ridgeloc.tolist()]),'%s/ridge.json' % (dbs_path))
    for bound in boundtypes:
        #name, segs = pickle.load(open('%s/%s.pkl' % (dbs_path,bound), 'rb'))
        
        try:
            name, segs = tools.readjson('%s/%s.json' % (dbs_path,bound))
        except IOError: #Download to default directory
            filedir = tools.get_filedir(checkwrite=True,makedir=True)
            data.update_file('%s.json' % (bound))
            name, segs = tools.readjson('%s/%s.json' % (filedir,bound))
        
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

<<<<<<< HEAD
def globalmap(ax,valarray,vmin,vmax,dbs_path='.',colorlabel=None,colorticks=True,colorpalette='rem3d',colorcontour=20,hotspots=False,grid=[30.,90.],gridwidth=0, **kwargs):
=======

def plot_gcpaths(m,stlon,stlat,eplon,eplat,ifglobal=True,**kwargs):
    """plots great-circle paths from lon lat arrays. Uses cartopy."""
    if kwargs:
        for x,y,z,w in np.vstack((stlon,stlat,eqlon,eqlat)).transpose():
            m.drawgreatcircle(x,y,z,w, **kwargs)
        x, y = m(stlon, stlat); m.scatter(x, y, marker='^', **kwargs)
        x, y = m(eplon, eplat); m.scatter(x, y, marker='o', **kwargs)
    else:
        for x,y,z,w in np.vstack((stlon,stlat,eqlon,eqlat)).transpose():
            m.drawgreatcircle(x,y,z,w)
        x, y = m(stlon, stlat); m.scatter(x, y, marker='^', edgecolors='k')
        x, y = m(eplon, eplat); m.scatter(x, y, marker='o', edgecolors='k')
    m.coastlines(color='gray')
    if ifglobal: m.set_global()    # set global extent
    return m


def backgroundmap(ax,dbs_path='.',platescolor='r', **kwargs):
    """plots a background map of a 3D model on axis ax. kwargs are arguments for Basemap"""
    
    # set up map
    if kwargs:
        m = Basemap(ax=ax, **kwargs)
    else:
        m = Basemap(ax=ax,projection='robin', lat_0=0, lon_0=150, resolution='l')
    
    clip_path = m.drawmapboundary()
    # draw coastlines.
#     m.drawcoastlines(linewidth=1.)
    # draw a boundary around the map, fill the background.
    # this background will end up being the ocean color, since
    # the continents will be drawn on top.
    m.drawmapboundary(fill_color='white')
    # fill continents, set lake color same as ocean color.
    m.fillcontinents(color='lightgray',lake_color='white')
    # add plates and hotspots
    dbs_path=tools.get_fullpath(dbs_path)
    plot_plates(m, dbs_path=dbs_path, color=platescolor, linewidth=1.)
    m.drawmapboundary(linewidth=1.5)    
    return m

def insetgcpathmap(ax,lat1,lon1,azimuth,gcdelta,projection='ortho',width=50.,height=50.,dbs_path='.',platescolor='r',numdegticks=7,hotspots=False):
    """plots the great-circle path between loc1-loc2. takes width/heght arguments in degrees if proj is merrcator,etc."""
    
    # Calculate intermediate points    
    lat2,lon2=mapping.getDestinationLatLong(lat1,lon1,azimuth,gcdelta*111325.)
    interval=gcdelta*111325./(numdegticks-1) # interval in km
    coords=np.array(mapping.getintermediateLatLong(lat1,lon1,azimuth,gcdelta*111325.,interval))

    # Center lat lon based on azimuth
    if gcdelta > 350.:
        lat_0,lon_0=mapping.getDestinationLatLong(lat1,lon1,azimuth,45.*111325.)
    elif gcdelta >= 180. and gcdelta <= 350.:
        lat_0,lon_0=mapping.getDestinationLatLong(lat1,lon1,azimuth,90.*111325.)
    else:
        lat_0,lon_0=mapping.getDestinationLatLong(lat1,lon1,azimuth,gcdelta/2.*111325.)
        
    # Choose what to do based on projection
    if projection=='ortho':
        m=backgroundmap(ax,tools.get_fullpath(dbs_path),projection=projection, lat_0=lat_0, lon_0=lon_0, resolution='l')
    else:
        # center left lat/lon, then left crnr
        latcenleft,loncenleft=mapping.getDestinationLatLong(lat_0,lon_0,-90.,width*111325./2.)
        llcrnrlat,llcrnrlon=mapping.getDestinationLatLong(latcenleft,loncenleft,180.,height*111325./2.)
        # center right lat/lon, then left crnr
        latcenright,loncenright=mapping.getDestinationLatLong(lat_0,lon_0,90.,width*111325./2.)
        urcrnrlat,urcrnrlon=mapping.getDestinationLatLong(latcenright,loncenright,0.,height*111325./2.)

        m=backgroundmap(ax,tools.get_fullpath(dbs_path),projection=projection, lat_0=lat_0, lon_0=lon_0, resolution='l',llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
        # draw parallels and meridians.
        # label parallels on right and top
        # meridians on bottom and left
        parallels = np.arange(-90,91,10.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels,labels=[1,0,0,0])
        meridians = np.arange(0.,361.,20.)
        m.drawmeridians(meridians,labels=[0,0,0,1])

    if hotspots: plot_hotspots(m, dbs_path=tools.get_fullpath(dbs_path), s=30, color='lightgreen', edgecolor='k')

    if gcdelta > 350.:
        dotsstart_x,dotsstart_y=m(coords[0:1][:,1],coords[0:1][:,0])
        m.scatter(dotsstart_x,dotsstart_y,s=50,zorder=10,facecolor='orange',edgecolor='k')
    dotsstart_x,dotsstart_y=m(coords[1:2][:,1],coords[1:2][:,0])
    dots_x,dots_y=m(coords[2:-1][:,1],coords[2:-1][:,0])
    m.scatter(dots_x,dots_y,s=50,zorder=10,facecolor='w',edgecolor='k')
    m.scatter(dotsstart_x,dotsstart_y,s=50,zorder=10,facecolor='m',edgecolor='k')
    dotsall_x,dotsall_y=m(coords[:,1],coords[:,0])
    if gcdelta < 180.:
        m.drawgreatcircle(lon1, lat1, lon2, lat2,color='k',linewidth=3.)
    elif gcdelta == 180.:
        latextent1,lonextent1=mapping.getDestinationLatLong(lat1,lon1,azimuth,1.*111325.)
        latextent2,lonextent2=mapping.getDestinationLatLong(lat1,lon1,azimuth,178.*111325.)
#         latextent2,lonextent2=mapping.getDestinationLatLong(lat_0,lon_0,180.+azimuth,89.*111325.)
        lonextent,latextent=m([lonextent1,lonextent2],[latextent1,latextent2])
        m.plot(lonextent,latextent,color='k',linewidth=3.)
    return m


def gettopotransect(lat1,lng1,azimuth,gcdelta,filename='ETOPO1_Bed_g_gmt4.grd',dbs_path='.',plottopo=False,numeval=50,downsampleetopo1=True,distnearthreshold=5.,outfile='transect_topo.json'):
    """Get the topography transect.dbs_path should have filename. numeval is number of evaluations of topo/bathymetry along the transect. downsampleetopo1 if true samples every 3 points to get 0.5 deg before interpolation. distnearthreshold is the threshold for nearest points to consider for interpolation in km. """
    
    recalculate=True
    if os.path.isfile(outfile):
        [lat1_t,lng1_t,azimuth_t,gcdelta_t,filename_t,dbs_path_t,plottopo_t,numeval_t,downsampleetopo1_t,distnearthreshold_t],evalpoints_t,grid_z1_t=tools.readjson(outfile)
        dbs_path=tools.get_fullpath(dbs_path)
        if lat1_t==lat1 and lng1_t==lng1 and azimuth_t==azimuth and gcdelta_t==gcdelta and str(filename_t)==filename and str(dbs_path_t)==dbs_path and plottopo_t==plottopo and numeval_t==numeval and downsampleetopo1_t==downsampleetopo1 and distnearthreshold_t==distnearthreshold:
            evalpoints=evalpoints_t
            grid_z1=grid_z1_t
            recalculate=False
    
    if recalculate:    
        # Calculate intermediate points            
        lat2,lng2=mapping.getDestinationLatLong(lat1,lng1,azimuth,gcdelta*111325.)
        interval=gcdelta*111325./(numeval-1) # interval in km
        coords=np.array(mapping.getintermediateLatLong(lat1,lng1,azimuth,gcdelta*111325.,interval))
        
        if(len(coords) != numeval):
            print "Error: The number of intermediate points is not accurate. Decrease it?"
            pdb.set_trace()
            sys.exit(2)
    
        # Read the topography file
        dbs_path=tools.get_fullpath(dbs_path)
        if os.path.isfile(dbs_path+'/'+filename):
            data = netcdf(dbs_path+'/'+filename,'r')
        else:
            print "Error: Could not find file "+dbs_path+'/'+filename    
            pdb.set_trace()
            sys.exit(2)
    
        # Create the arrays for interpolation
        radarr=np.array([6371.])
        if downsampleetopo1: 
            sampling=3
        else:
            sampling=1    
        lonarr=data.variables['lon'][::][0::sampling] #sample every 3 points to get 0.5 deg
        latarr=data.variables['lat'][::][0::sampling] #sample every 3 points to get 0.5 deg
        grid_y, grid_z, grid_x=np.meshgrid(latarr,radarr,lonarr)
    #     grid_x, grid_y = np.meshgrid(data.variables['lon'][::],data.variables['lat'][::])    
        values=data.variables['z'][::][0::sampling,0::sampling] #sample every 3 points to get 0.5 deg    
        if plottopo:
            fig=plt.figure() 
            ax=fig.add_subplot(1,1,1)
            ax.pcolormesh(grid_x,grid_y,values,cmap='coolwarm')
            plt.show()
        #Evaluate the topography at the points
        rlatlon=np.column_stack((np.ravel(grid_z), np.ravel(grid_y), np.ravel(grid_x)))

        t0 = time.time()
        gridpoints=mapping.spher2cart(rlatlon)
        evalpoints=np.column_stack((6371.*np.ones_like(coords[:,1]),coords[:,0],coords[:,1]))
        coordstack=mapping.spher2cart(evalpoints)
        checkifnear=np.zeros_like(gridpoints[:,0],dtype=bool) # array for checking if the point is near to any point in the path (coordstack)
    
        print "....Getting the topography transect from "+dbs_path+'/'+filename
        bar = progressbar.ProgressBar(maxval=len(coordstack), \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        for ii in np.arange(len(coordstack)):
            dist=np.sqrt( np.power(gridpoints[:,0]-coordstack[ii,0],2) + np.power(gridpoints[:,1]-coordstack[ii,1],2) + np.power(gridpoints[:,2]-coordstack[ii,2],2) )
            checkifneartemp=dist<distnearthreshold
            checkifnear=np.logical_or(checkifnear,checkifneartemp)
            if len(np.where(checkifneartemp==True)[0]) == 0: 
                print "Error: No points found with the given distnearthreshold in atleast 1 point. Increase the bound? "
                sys.exit(2)
            bar.update(ii+1)    
        bar.finish()
    
        indexselect=np.where(checkifnear==True)[0]    
        valuesselect=np.ravel(values)[indexselect]    
        print("--- Chose %s points with distances < %s km for interpolation ---" % (len(valuesselect),  distnearthreshold))
        grid_z1 = spint.griddata(gridpoints[indexselect], valuesselect, coordstack, method='nearest')
    
        writearr=np.array([[lat1,lng1,azimuth,gcdelta,filename,dbs_path,plottopo,numeval,downsampleetopo1,distnearthreshold],evalpoints.tolist(),grid_z1.tolist()])
        tools.writejson(writearr,outfile)
        data.close() #close netcdf file
    return evalpoints,grid_z1

def globalmap(ax,valarray,vmin,vmax,dbs_path='.',colorlabel=None,colorpalette='rem3d',colorcontour=21,hotspots=False,grid=[30.,90.],gridwidth=0, **kwargs):
>>>>>>> 8a2055a34aa4bcbebba0fd840df18bd5cb5cde75
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
    
    colorticks : Label and draw the ticks in the colorbar if True (default)
    
    projection : map projection for the global plot
    
    colorlabel : label to use for the colorbar. If None, no colorbar is plotted.
    
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
    m.drawparallels(parallels,labels=[True,False,False,False],linewidth=gridwidth)
    meridians = np.arange(-180.,180.,grid[1])
    m.drawmeridians(meridians,labels=[False,False,False,True],linewidth=gridwidth)

    # Get the color map
    try:
        cpalette = plt.get_cmap(colorpalette)
    except ValueError:
        cpalette=standardcolorpalette(RGBoption=colorpalette)
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
    elif isinstance(colorcontour,(int, float)): # Number of intervals for color bar
        bounds = np.linspace(vmin,vmax,colorcontour+1)
        mytks = np.arange(vmin,vmax+(vmax-vmin)/4.,(vmax-vmin)/4.)
        mytkslabel = [str(a) for a in mytks]
        spacing='proportional'
    else:
        print("Error: Undefined colorcontour in globalmap; should be a numpy array, list or integer ")
        sys.exit(2)        
    norm = mcolors.BoundaryNorm(bounds,cpalette.N)
    
    # plot the model
    for ii in np.arange(len(valarray['lon'])): 
        if valarray['lon'][ii] > 180.: valarray['lon'][ii]=valarray['lon'][ii]-360.
    numlon=len(np.unique(valarray['lon']))
    numlat=len(np.unique(valarray['lat']))
    # grid spacing assuming a even grid
    # Get the unique lat and lon spacing, avoiding repeated lat/lon
    spacing_lat = np.ediff1d(np.sort(valarray['lat']))
    spacing_lat = np.unique(spacing_lat[spacing_lat != 0])
    spacing_lon = np.ediff1d(np.sort(valarray['lon']))
    spacing_lon = np.unique(spacing_lon[spacing_lon != 0])
    # Check if an unique grid spacing exists for both lat and lon
    if len(spacing_lon)!=1 or len(spacing_lat)!=1 or np.any(spacing_lat!=spacing_lon): 
        print("Warning: spacing for latitude and longitude should be the same. Using nearest neighbor interpolation")
        # compute native map projection coordinates of lat/lon grid.
        x, y = m(valarray['lon'], valarray['lat'])
        rlatlon = np.vstack([np.ones(len(valarray['lon'])),valarray['lat'],valarray['lon']]).transpose()
        xyz = mapping.spher2cart(rlatlon)
        
        # Create a grid
        grid_spacing = 1.
        lat = np.arange(-90.+grid_spacing/2.,90.+grid_spacing/2.,grid_spacing)
        lon = np.arange(-180.+grid_spacing/2.,180.+grid_spacing/2.,grid_spacing)
        lons,lats=np.meshgrid(lon,lat)
        
        # evaluate in cartesian
        rlatlon = np.vstack([np.ones_like(lons.flatten()),lats.flatten(),lons.flatten()]).transpose()
        xyz_new = mapping.spher2cart(rlatlon)
        
        # grid the data.
        val = spint.griddata(xyz, valarray['val'], xyz_new, method='nearest').reshape(lons.shape)    
        s = m.transform_scalar(val,lon,lat, 1000, 500)
        im = m.imshow(s, cmap=cpalette.name, vmin=vmin, vmax=vmax, norm=norm)
        #im = m.contourf(lons, lats,val, norm=norm, cmap=cpalette.name, vmin=vmin, vmax=vmax,latlon=True)
    else: 
        grid_spacing = spacing_lat
        # Create a grid
        lat = np.arange(-90.+grid_spacing/2.,90.+grid_spacing/2.,grid_spacing)
        lon = np.arange(-180.+grid_spacing/2.,180.+grid_spacing/2.,grid_spacing)
        X,Y=np.meshgrid(lon,lat)
        val = np.empty_like(X)
        val[:] = np.nan;
        for i in range(0, valarray['lat'].size):
            ilon = np.where(X[0,:]==valarray['lon'][i])[0][0]
            ilat = np.where(Y[:,0]==valarray['lat'][i])[0][0]
            val[ilat,ilon] = valarray['val'][i]
        s = m.transform_scalar(val,lon,lat, 1000, 500)
        im = m.imshow(s, cmap=cpalette.name, vmin=vmin, vmax=vmax, norm=norm)
    # add plates and hotspots
    dbs_path=tools.get_fullpath(dbs_path)
    plot_plates(m, dbs_path=dbs_path, color='w', linewidth=1.5)

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
        if colorticks:
            cbar.ax.tick_params(labelsize=14)
            # Customize colorbar tick labels
            cbar.set_ticks(mytks)
            mytkslabels = [str(int(a)) if (a).is_integer() else str(a) for a in mytks]
            cbar.ax.set_yticklabels(mytkslabels)
        else:   
            # Remove color bar tick lines, while keeping the tick labels
            cbarytks = plt.getp(cbar.ax.axes, 'yticklines')
            plt.setp(cbarytks, visible=False)
            cbarytks = plt.getp(cbar.ax.axes, 'yticklabels')
            plt.setp(cbarytks, visible=False)
        # Colorbar label, customize fontsize and distance to colorbar
        cbar.set_label(colorlabel,rotation=90, labelpad=5)

        
    m.drawmapboundary(linewidth=1.5)    
    if hotspots: plot_hotspots(m, dbs_path=dbs_path, s=30, color='m', edgecolor='k')

    return m    

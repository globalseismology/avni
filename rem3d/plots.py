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
from mpl_toolkits.basemap import Basemap
import multiprocessing
import cartopy.crs as ccrs
import codecs,json #printing output
from joblib import Parallel, delayed
import pdb	#for the debugger pdb.set_trace()
# from scipy.io import netcdf_file as netcdf #reading netcdf files
from netCDF4 import Dataset as netcdf #reading netcdf files
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import itertools
import time
import progressbar
# For polar sectionplot
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator,DictFormatter,FixedLocator
from matplotlib import gridspec # Relative size of subplots
from numba import jit

####################       IMPORT OWN MODULES     ######################################
from . import geolib
from . import tools
########################  	GENERIC   ################################################   	    		    
	    		    
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
    
def customcolorpalette(name='bk',cptfolder='~/CPT',colorlist=None,colormax=2.,middlelimit=0.5,ifgraytest=0):
	"""Used to return preset color palettes from cptfolder. ifgraytest test how the figure looks in gray scale. (-colormax,colormax) are the limits of the colorbar. zerolimit is the limit to which the middle color (e.g. grey) will extend on either side of colorttmax mid. """
	c = mcolors.ColorConverter().to_rgb	
	if name=='r_lgrey_b':
		colorlist=[c('blue'), c('lightgray'), (2.*colormax-2.*middlelimit)/(4.*colormax), c('lightgray'),c('lightgray'), (2.*colormax+2.*middlelimit)/(4.*colormax), c('lightgray'),c('red'), 1., c('red')]
	elif name=='bk':
		cptfolder=tools.get_fullpath(cptfolder)
		if os.path.isfile(cptfolder+'/bk1_0.cpt_'):
			colorlist=getcolorlist(cptfolder+'/bk1_0.cpt_')
		else:
			print "Error: Could not find file "+cptfolder+'/bk1_0.cpt_'	
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
	register_cmap(name=custom_cmap.name, cmap=custom_cmap)
	palette=custom_cmap.name
	
	if ifgraytest==1:
		palette=grayify_cmap(palette)
		
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


def plot_gcpaths(m,stlon,stlat,eplon,eplat,ifglobal=True,**kwargs):
	"""plots great-circle paths from lon lat arrays. Uses cartopy."""
	if kwargs:
		m.plot([stlon,eplon],[stlat,eplat], transform=ccrs.Geodetic(), **kwargs)
		m.scatter(stlon, stlat, transform=ccrs.Geodetic(), marker='^',edgecolors='k', **kwargs)
		m.scatter(eplon, eplat, transform=ccrs.Geodetic(), marker='o',edgecolors='k', **kwargs)
	else:
		m.plot([lon1,lon2],[lat1,lat2], transform=ccrs.Geodetic())
		m.scatter(stlon, stlat, transform=ccrs.Geodetic(), marker='^', edgecolors='k')
		m.scatter(eplon, eplat, transform=ccrs.Geodetic(), marker='o', edgecolors='k')
	m.coastlines(color='gray')
	if ifglobal: m.set_global()	# set global extent
	return m


def backgroundmap(ax,dbs_path='.',platescolor='r', **kwargs):
	"""plots a background map of a 3D model on axis ax. kwargs are arguments for Basemap"""
	
	# set up map
	if kwargs:
		m = Basemap(ax=ax, **kwargs)
	else:
		m = Basemap(ax=ax,projection='hammer', lat_0=0, lon_0=150, resolution='l')
	
	clip_path = m.drawmapboundary()
	# draw coastlines.
# 	m.drawcoastlines(linewidth=1.)
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
	lat2,lon2=geolib.getDestinationLatLong(lat1,lon1,azimuth,gcdelta*111325.)
	interval=gcdelta*111325./(numdegticks-1) # interval in km
	coords=np.array(geolib.getintermediateLatLong(lat1,lon1,azimuth,gcdelta*111325.,interval))

	# Center lat lon based on azimuth
	if gcdelta > 350.:
		lat_0,lon_0=geolib.getDestinationLatLong(lat1,lon1,azimuth,45.*111325.)
	elif gcdelta >= 180. and gcdelta <= 350.:
		lat_0,lon_0=geolib.getDestinationLatLong(lat1,lon1,azimuth,90.*111325.)
	else:
		lat_0,lon_0=geolib.getDestinationLatLong(lat1,lon1,azimuth,gcdelta/2.*111325.)
		
	# Choose what to do based on projection
	if projection=='ortho':
		m=backgroundmap(ax,tools.get_fullpath(dbs_path),projection=projection, lat_0=lat_0, lon_0=lon_0, resolution='l')
	else:
		# center left lat/lon, then left crnr
		latcenleft,loncenleft=geolib.getDestinationLatLong(lat_0,lon_0,-90.,width*111325./2.)
		llcrnrlat,llcrnrlon=geolib.getDestinationLatLong(latcenleft,loncenleft,180.,height*111325./2.)
		# center right lat/lon, then left crnr
		latcenright,loncenright=geolib.getDestinationLatLong(lat_0,lon_0,90.,width*111325./2.)
		urcrnrlat,urcrnrlon=geolib.getDestinationLatLong(latcenright,loncenright,0.,height*111325./2.)

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
		latextent1,lonextent1=geolib.getDestinationLatLong(lat1,lon1,azimuth,1.*111325.)
		latextent2,lonextent2=geolib.getDestinationLatLong(lat1,lon1,azimuth,178.*111325.)
# 		latextent2,lonextent2=geolib.getDestinationLatLong(lat_0,lon_0,180.+azimuth,89.*111325.)
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
		lat2,lng2=geolib.getDestinationLatLong(lat1,lng1,azimuth,gcdelta*111325.)
		interval=gcdelta*111325./(numeval-1) # interval in km
		coords=np.array(geolib.getintermediateLatLong(lat1,lng1,azimuth,gcdelta*111325.,interval))
		
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
	# 	grid_x, grid_y = np.meshgrid(data.variables['lon'][::],data.variables['lat'][::])	
		values=data.variables['z'][::][0::sampling,0::sampling] #sample every 3 points to get 0.5 deg	
		if plottopo:
			fig=plt.figure() 
			ax=fig.add_subplot(1,1,1)
			ax.pcolormesh(grid_x,grid_y,values,cmap='coolwarm')
			plt.show()
		#Evaluate the topography at the points
		rlatlon=np.column_stack((np.ravel(grid_z), np.ravel(grid_y), np.ravel(grid_x)))

		t0 = time.time()
		gridpoints=geolib.spher2cart(rlatlon)
		evalpoints=np.column_stack((6371.*np.ones_like(coords[:,1]),coords[:,0],coords[:,1]))
		coordstack=geolib.spher2cart(evalpoints)
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

def globalmap(ax,valarray,vmin,vmax,dbs_path='.',colorlabel=None,colorpalette='bk',colorcontour=20,hotspots=False,grid=[30.,90.],gridwidth=0, **kwargs):
	"""plots a 2-D cross-section of a 3D model on axis ax. kwargs are arguments for Basemap. color* are for the colormap used."""

	# set up map
	if kwargs:
		m = Basemap(ax=ax, **kwargs)
	else:
		m = Basemap(ax=ax,projection='hammer', lat_0=0, lon_0=150, resolution='c')
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
	lontemp = np.transpose(valarray['lon'].reshape((numlon,numlat)))
	lonsort = np.argsort(lontemp) # Sort the longitude for transform_scalar
	lon = np.sort(lontemp) # Sort the longitude for transform_scalar
	lattemp = np.transpose(valarray['lat'].reshape((numlon,numlat)))
	valtemp = np.transpose(valarray['val'].reshape((numlon,numlat)))
	val=np.zeros_like(valtemp);lat=np.zeros_like(lattemp)
	for ii in np.arange(numlat):
		for jj in np.arange(numlon):
			val[ii,jj]=valtemp[ii,lonsort[ii,jj]]
			lat[ii,jj]=lattemp[ii,lonsort[ii,jj]]
	s = m.transform_scalar(val,lon[0,:], lat[:,0], 1000, 500)
	# Get the color map
# 	c = mcolors.ColorConverter().to_rgb
# 	colorlist=[c('red'), c('violet'), 0.33, c('violet'), c('blue'), 0.66, c('blue')]
# 	cpalette=customcolorpalette('custom',colorlist)
	cpalette=customcolorpalette(colorpalette)
	# define the 10 bins and normalize
	bounds = np.linspace(vmin,vmax,colorcontour+1)
	norm = mcolors.BoundaryNorm(bounds,cpalette.N)
	im = m.imshow(s, cmap=cpalette.name, clip_path=clip_path, vmin=vmin, vmax=vmax, norm=norm)
#	# add plates and hotspots
	dbs_path=tools.get_fullpath(dbs_path)
	plot_plates(m, dbs_path=dbs_path, color='w', linewidth=1.5)
	if hotspots: plot_hotspots(m, dbs_path=dbs_path, s=30, color='m', edgecolor='k')

# add a colorbar
	if colorlabel is not None:
# 		cb = plt.colorbar(im,orientation='vertical',fraction=0.05,pad=0.05)
# 		cb.set_label(colorlabel)
		# Set colorbar, aspect ratio
		cbar = plt.colorbar(im, alpha=0.05, aspect=12, shrink=0.4,norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds,extendrect=True)
		cbar.solids.set_edgecolor("face")
		# Remove colorbar container frame
# 		cbar.outline.set_visible(False)
		# Fontsize for colorbar ticklabels
		cbar.ax.tick_params(labelsize=14)
		# Customize colorbar tick labels
		mytks = np.arange(vmin,vmax+(vmax-vmin)/4.,(vmax-vmin)/4.)
		cbar.set_ticks(mytks)
		cbar.ax.set_yticklabels([str(a) for a in mytks])
		# Colorbar label, customize fontsize and distance to colorbar
		cbar.set_label(colorlabel,rotation=90, fontsize=16, labelpad=5)
		# Remove color bar tick lines, while keeping the tick labels
# 		cbarytks = plt.getp(cbar.ax.axes, 'yticklines')
# 		plt.setp(cbarytks, visible=False)
		
	m.drawmapboundary(linewidth=1.5)	
	return m

def setup_axes(fig, rect, theta, radius, numdegticks=7,r_locs = [3480.,3871.,4371.,4871.,5371.,5871.,6346.6],r_labels = ['CMB',' ','2000',' ','1000',' ','Moho'],fontsize=12):
	"""Setup the polar axis for section plot. numdegticks is the number of grids in theta. rect is the 3-digit number for axis on a plot. (theta, radius) are array for the range."""
    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
	tr = Affine2D().scale(np.pi/180., 1.) + PolarAxes.PolarTransform()

    # Find grid values appropriate for the coordinate (degree).
    # The argument is an approximate number of grids.
#     theta_grid_locator = angle_helper.LocatorD(numdegticks)

	theta_grid_locator=FixedLocator(np.arange(theta[0], theta[1], numdegticks))
    # And also use an appropriate formatter:
	theta_tick_formatter = angle_helper.FormatterDMS()

    # set up number of ticks for the r-axis
#     r_grid_locator = MaxNLocator(7)
	r_grid_locator=FixedLocator(r_locs)
	
    # Plot the radius ticks    
	r_ticks = {loc : label for loc, label in zip(r_locs, r_labels)}
	r_tick_formatter = DictFormatter(r_ticks)
	
	# the extremes are passed to the function
	grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                extremes=(theta[0], theta[1], radius[0], radius[1]),
                                grid_locator1=theta_grid_locator,
                                grid_locator2=r_grid_locator,
                                tick_formatter1=theta_tick_formatter,
                                tick_formatter2=r_tick_formatter,
                                )
	
	ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
	fig.add_subplot(ax1)
	
	if theta[1]-theta[0]>350.:
		ax1.axis["bottom"].set_visible(False)
		ax1.axis["top"].set_visible(False)
		ax1.axis["left"].set_visible(False)
		ax1.axis["right"].set_visible(False)
	elif theta[1]-theta[0]>90. and theta[1]-theta[0]<=180. :
		# Make ticks on outside
		ax1.axis["left"].set_axis_direction("top")
		ax1.axis["right"].set_axis_direction("bottom")
		# Make ticks invisible on top and bottom axes
		ax1.axis["bottom"].set_visible(False)
		ax1.axis["top"].set_visible(False)
		ax1.axis["left"].toggle(ticklabels=False, label=False)
		ax1.axis["right"].toggle(ticklabels=False, label=False)
	else:
		# adjust axis
		# the axis artist lets you call axis with
		# "bottom", "top", "left", "right"
		# Make ticks on outside
		ax1.axis["left"].set_axis_direction("top")
		ax1.axis["right"].set_axis_direction("top")
		
		# Make ticks invisible on top and bottom axes
		ax1.axis["bottom"].set_visible(False)
		ax1.axis["top"].set_visible(False)
	#     ax1.axis["top"].set_axis_direction("bottom")
	#     ax1.axis["top"].toggle(ticklabels=False, label=False)
	#     ax1.axis["top"].major_ticklabels.set_axis_direction("top")
	#     ax1.axis["top"].label.set_axis_direction("top")
	# 
		ax1.axis["left"].toggle(ticklabels=False, label=False)
		ax1.axis["right"].toggle(ticklabels=True, label=True)
		ax1.axis["right"].label.set_axis_direction("bottom")
		ax1.axis["right"].major_ticklabels.set_axis_direction("bottom")
		ax1.axis["right"].major_ticklabels.set_fontsize(fontsize)
		ax1.axis["right"].label.set_text("Depth (km)")
		ax1.axis["right"].label.set_fontsize(fontsize)
	#     ax1.axis["top"].label.set_text(ur"$\alpha$ [\u00b0]")

	# create a parasite axes
	aux_ax = ax1.get_aux_axes(tr)
	
	aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
	ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                         # drawn twice, and possibly over some other
                         # artists. So, we decrease the zorder a bit to
                         # prevent this.

	# plot the degree increments at 6371 km, always on top (large zorder)
	degticks=np.linspace(theta[0],theta[1],numdegticks)
	if theta[1]-theta[0]>350.:
		degticks0=degticks[0:1] #color first tick in magenta
		aux_ax.scatter(degticks0,6346.6*np.ones(len(degticks0)),s=50,clip_on=False,zorder=10,facecolor='orange',edgecolor='k')

		degticksstart=degticks[1:2] #color first tick in magenta
		degticks=degticks[2:-1] # do not not plot the first and last 2 ticks
	else:
		degticksstart=degticks[1:2] #color first tick in magenta
		degticks=degticks[2:-1] # do not not plot the first and last 2 ticks
	aux_ax.scatter(degticks,6346.6*np.ones(len(degticks)),s=50,clip_on=False,zorder=10,facecolor='w',edgecolor='k')
	aux_ax.scatter(degticksstart,6346.6*np.ones(len(degticksstart)),s=50,clip_on=False,zorder=10,facecolor='m',edgecolor='k')
    
	# 410 and 650 
	theta_arr = np.linspace(theta[0],theta[1])
	disc_arr=[5961.,5721.]
	for disc in disc_arr:
		aux_ax.plot(theta_arr, disc*np.ones(len(theta_arr)),'k--', linewidth=1.5,zorder=10)
	# Surface ICB CMB
	# 	disc_arr=[6371.,3480.,1215.]
	disc_arr=[6346.6,3480.]
	for disc in disc_arr:
		aux_ax.plot(theta_arr, disc*np.ones(len(theta_arr)), 'k', linewidth=1,zorder=9)
	
	return ax1, aux_ax
    
def plottopotransect(ax,theta_range,elev,vexaggerate=150):
	"""Plot a section on the axis ax. """		
	elevplot1=np.array(elev)
	elevplot2=np.array(elev)
	# Blue for areas below sea level
	if np.min(elev)<0.:
		lowerlimit=6371.-np.min(elev)/1000.*vexaggerate
		elevplot2[elevplot2>0.]=0.
		ax.fill_between(theta_range, lowerlimit*np.ones(len(theta_range)),lowerlimit*np.ones(len(theta_range))+elevplot2/1000.*vexaggerate,facecolor='aqua', alpha=0.5)
	else:
		lowerlimit=6371.

	# Grey for areas above sea level
	elevplot1[elevplot1<0.]=0.
	ax.fill_between(theta_range, lowerlimit*np.ones(len(theta_range)), lowerlimit*np.ones(len(theta_range))+elevplot1/1000.*vexaggerate, facecolor='grey', alpha=0.5)

# 	title(phase, fontsize=20,loc='left')
	return ax
	

def getmodeltransect(lat1,lng1,azimuth,gcdelta,parameter='vs',radii=[3480.,6346.6],filename='S362ANI+M_kmps.nc',dbs_path='.',numeval=200,distnearthreshold=1000.,outfile='transect_model3D.json'):
	"""Get the topography transect.dbs_path should have filename. numeval is number of evaluations of topo/bathymetry along the transect. downsampleetopo1 if true samples every 3 points to get 0.5 deg before interpolation. distnearthreshold is the threshold for nearest points to consider for interpolation in km. """
	
	recalculate=True
	if os.path.isfile(outfile):
		[lat1_t,lng1_t,azimuth_t,gcdelta_t,parameter_t,filename_t,dbs_path_t,numeval_t,distnearthreshold_t],radii_t,evalpoints_t,grid_z1_t=tools.readjson(outfile)
		dbs_path=tools.get_fullpath(dbs_path)
		if lat1_t==lat1 and lng1_t==lng1 and azimuth_t==azimuth and gcdelta_t==gcdelta and str(parameter_t)==parameter and str(filename_t)==filename and str(dbs_path_t)==dbs_path and numeval_t==numeval and distnearthreshold_t==distnearthreshold and radii_t ==radii:
			evalpoints=evalpoints_t
			grid_z1=grid_z1_t
			recalculate=False
	
	if recalculate:	
		# Calculate intermediate points			
		lat2,lng2=geolib.getDestinationLatLong(lat1,lng1,azimuth,gcdelta*111325.)
		interval=gcdelta*111325./(numeval-1) # interval in km
		radevalarr=np.linspace(radii[0],radii[1],numeval) #radius arr in km
		coords=np.array(geolib.getintermediateLatLong(lat1,lng1,azimuth,gcdelta*111325.,interval))
		
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
		radarr=6371.-data.variables['depth'][::]	
		lonarr=data.variables['longitude'][::]
		latarr=data.variables['latitude'][::]
		grid_y, grid_z,grid_x=np.meshgrid(latarr,radarr,lonarr)
		values=data.variables[parameter][::]
		#Evaluate the topography at the points
		rlatlon=np.column_stack((np.ravel(grid_z), np.ravel(grid_y), np.ravel(grid_x)))

		t0 = time.time()
		gridpoints=geolib.spher2cart(rlatlon)
		evalpoints=np.column_stack((radevalarr[0]*np.ones_like(coords[:,1]),coords[:,0],coords[:,1]))
		for radius in radevalarr[1:]:
			pointstemp=np.column_stack((radius*np.ones_like(coords[:,1]),coords[:,0],coords[:,1]))
			evalpoints=np.row_stack((evalpoints,pointstemp))	

		coordstack=geolib.spher2cart(evalpoints)
		checkifnear=np.zeros_like(gridpoints[:,0],dtype=bool) # array for checking if the point is near to any point in the path (coordstack)
	
		
		print "....Getting the 3D-model transect from "+dbs_path+'/'+filename
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
	
		writearr=np.array([[lat1,lng1,azimuth,gcdelta,parameter,filename,dbs_path,numeval,distnearthreshold],radii,evalpoints.tolist(),grid_z1.tolist()])
				
		tools.writejson(writearr,outfile)
		data.close() #close netcdf file
	return evalpoints,grid_z1
	

def plot1section(lat1,lng1,azimuth,gcdelta,model=None,vmin=None,vmax=None,dbs_path='.', colorlabel=None,colorpalette='bk',colorcontour=20,nelevinter=50,radii=[3480.,6346.6],n3dmodelinter=200,vexaggerate=150,figuresize=[8,4],width_ratios=[1, 2],outfile=None):
	"""Plot one section through the Earth through a pair of points.""" 
	
	fig = plt.figure(1, figsize=(figuresize[0],figuresize[1]))
	gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios) 

	fig.subplots_adjust(wspace=0.3, left=0.05, right=0.95)
	fig.patch.set_facecolor('white')
	
	# Specify theta such that it is symmetric
	lat2,lng2=geolib.getDestinationLatLong(lat1,lng1,azimuth,gcdelta*111325.)
	if gcdelta==180. or gcdelta==360.:
		theta=[0.,gcdelta]
	else:
		theta=[90.-gcdelta/2.,90.+gcdelta/2.]
	theta_range=np.linspace(theta[0],theta[1],nelevinter)
	# default is not to extend radius unless vexaggerate!=0
	extend_radius=0.
	if vexaggerate != 0:
		rlatlon,elev=gettopotransect(lat1,lng1,azimuth,gcdelta,dbs_path=dbs_path,numeval=nelevinter)
		if min(elev)< 0.:
	 		extend_radius=(max(elev)-min(elev))*vexaggerate/1000.
		else:
	 		extend_radius=max(elev)*vexaggerate/1000.
	
	####### Inset map
	ax=plt.subplot(gs[0])
	if gcdelta > 180.:
		numdegticks=13
		insetgcpathmap(ax,lat1,lng1,azimuth,gcdelta,projection='ortho',dbs_path=dbs_path,numdegticks=numdegticks)
	elif gcdelta >= 30. and gcdelta <=180:
		numdegticks=7
		insetgcpathmap(ax,lat1,lng1,azimuth,gcdelta,projection='ortho',dbs_path=dbs_path,		numdegticks=numdegticks)
	else:	
		numdegticks=7
		width=gcdelta*1.2
		height=gcdelta*1.2
		insetgcpathmap(ax,lat1,lng1,azimuth,gcdelta,projection='merc',dbs_path=dbs_path,width=width,height=height,numdegticks=numdegticks)
	###### Actual cross-section
	ax1, aux_ax1 = setup_axes(fig, gs[1], theta, radius=[3480., 6371.+extend_radius],numdegticks=numdegticks)
	if vexaggerate != 0:
		aux_ax1=plottopotransect(aux_ax1,theta_range,elev,vexaggerate=vexaggerate)

	# Plot the model section
# 	grid_x, grid_y = np.mgrid[theta[0]:theta[1]:200j,3480.:6346.6:200j]
	grid_x, grid_y = np.meshgrid(np.linspace(theta[0],theta[1],n3dmodelinter),np.linspace(radii[0],radii[1],n3dmodelinter))
	junk,values = getmodeltransect(lat1,lng1,azimuth,gcdelta,radii=radii,dbs_path=dbs_path,numeval=n3dmodelinter)
	interp_values=np.array(values)
	
	interp_values=(interp_values-interp_values.mean()).reshape((n3dmodelinter,n3dmodelinter))
	#pdb.set_trace()
# 	values = np.random.rand(len(grid_x), len(grid_y))
	cpalette=customcolorpalette(colorpalette)
	# define the 10 bins and normalize
	bounds = np.linspace(vmin,vmax,colorcontour+1)
	norm = mcolors.BoundaryNorm(bounds,cpalette.N)
# 	im = m.imshow(s, cmap=cpalette.name, vmin=vmin, vmax=vmax, norm=norm)
	im=aux_ax1.pcolormesh(grid_x,grid_y,interp_values,cmap=cpalette.name,vmin=vmin, vmax=vmax, norm=norm)
	# add a colorbar
	if colorlabel is not None:
# 		cb = plt.colorbar(im,orientation='vertical',fraction=0.05,pad=0.05)
# 		cb.set_label(colorlabel)
		# Set colorbar, aspect ratio
		cbar = plt.colorbar(im, alpha=0.05, aspect=12, shrink=0.4,norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds,extendrect=True)
		cbar.solids.set_edgecolor('face')
		# Remove colorbar container frame
# 		cbar.outline.set_visible(False)
		# Fontsize for colorbar ticklabels
		cbar.ax.tick_params(labelsize=12)
		# Customize colorbar tick labels
		mytks = np.arange(vmin,vmax+(vmax-vmin)/4.,(vmax-vmin)/4.)
		cbar.set_ticks(mytks)
		cbar.ax.set_yticklabels([str(a) for a in mytks])
		# Colorbar label, customize fontsize and distance to colorbar
		cbar.set_label(colorlabel,rotation=90, fontsize=14, labelpad=5)
		# Remove color bar tick lines, while keeping the tick labels
# 		cbarytks = plt.getp(cbar.ax.axes, 'yticklines')
# 		plt.setp(cbarytks, visible=False)

	fig1 = plt.gcf()
	plt.show()
	plt.draw()
	if outfile is not None:
		fig1.savefig(outfile,dpi=100)
	return 

def plot1globalmap(filename,vmin,vmax,dbs_path='.',projection='hammer',colorlabel=None,lat_0=0,lon_0=150):
	"""Plot one global map""" 
	fig=plt.figure() 
	ax=fig.add_subplot(1,1,1)
	epixarr=models.readepixfile(filename)
	if 	projection=='ortho':
		globalmap(ax,epixarr,vmin,vmax,dbs_path,colorlabel,grid=[30.,30.],gridwidth=1,projection=projection,lat_0=lat_0, lon_0=lon_0)
	else:
		globalmap(ax,epixarr,vmin,vmax,dbs_path,colorlabel,grid=[30.,90.],gridwidth=0,projection=projection,lat_0=lat_0, lon_0=lon_0)
	plt.show()
	return 
	
	
	
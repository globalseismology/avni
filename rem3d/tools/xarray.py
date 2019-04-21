#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import os
import numpy as np
import xarray as xr
from scipy.spatial import cKDTree
import pickle

####################### IMPORT REM3D LIBRARIES  #######################################
from .trigd import sind
from ..mapping import spher2cart
from .. import constants
from .common import precision_and_scale,convert2nparray
#######################################################################################

def tree3D(treefile,latitude,longitude,radius_in_km):
    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)
    radius_in_km = convert2nparray(radius_in_km)

    #Build the tree if none is provided
    if os.path.isfile(treefile):
        print('... Reading KDtree file '+treefile)
        tree = pickle.load(open(treefile,'rb'))
    else:
        print('... KDtree file '+treefile+' not found for interpolations. Building it')
        rlatlon = np.column_stack((radius_in_km.flatten(),latitude.flatten(), longitude.flatten()))
        xyz = spher2cart(rlatlon)
        tree = cKDTree(xyz)
        pickle.dump(tree,open(treefile,'wb'))
    return tree

def querytree3D(tree,latitude,longitude,radius_in_km,values,nearest=1):    
    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)
    radius_in_km = convert2nparray(radius_in_km)
    evalpoints = np.column_stack((radius_in_km.flatten(),latitude.flatten(), longitude.flatten()))
    coordstack = spher2cart(evalpoints)
    d,inds = tree.query(coordstack,k=nearest)
    if nearest == 1:
        interp = values[inds]
    else:
        w = 1.0 / d**2
        interp = np.sum(w * values[inds], axis = 1)/ np.sum(w, axis=1)
    return interp 
    
def ncfile2tree3D(ncfile,treefile,lonlatdepth = ['longitude','latitude','depth'],stride=None, radius_in_km = None):
    """
    Read or write a pickle interpolant with KDTree

    Input Parameters:
    ----------------

    ncfile: netcdf format file e.g. topography
    
    stride: is the downsampling before interpolation. 
    
    radius_in_km: radius in kilometer when a 2D surface is valid. Ignores the 
                    3rd field in lonlatdepth.Typically 6371km for Earth

    lonlatdepth: variable name of the longitude, latitude, depth (in km) arrays
    """

    #read topography file
    if os.path.isfile(ncfile):
        f = xr.open_dataset(ncfile)
    else:
        raise ValueError("Error: Could not find file "+ncfile)
    if stride != None:
        lon = f.variables[lonlatdepth[0]][::stride]
        lat = f.variables[lonlatdepth[1]][::stride]
    else:
        lon = f.variables[lonlatdepth[0]]
        lat = f.variables[lonlatdepth[1]]   
    if radius_in_km == None:
        if stride != None:
            dep = f.variables[lonlatdepth[2]][::stride]
        else:
            dep = f.variables[lonlatdepth[2]]
        rad = constants.R/1000. - dep
    else:
        rad = xr.IndexVariable('rad',[radius_in_km])
    f.close() #close netcdf file
    
    # get the tree
    gridlat, gridrad, gridlon = np.meshgrid(lat.data,rad.data,lon.data)
    tree = tree3D(treefile,gridlat,gridlon,gridrad)
    return tree


def AreaDataArray(data,latname = 'latitude', lonname = 'longitude'):
    """
    weighted average for xray data geographically averaged

    Parameters
    ----------
    data : Dataset or DataArray
        the xray object to average over

    Returns
    -------
    
    area: a DataArray object with area of each pixel
    
    """

    if not isinstance(data, xr.DataArray): raise ValueError("date must be an xray DataArray")
    
    pix_lat = np.unique(np.ediff1d(np.sort(data.coords[latname].values)))
    pix_lon = np.unique(np.ediff1d(np.sort(data.coords[lonname].values)))
    assert(len(pix_lat)==len(pix_lon)==1),'only one pixel size allowed in xarray'
    assert(pix_lat.item()==pix_lon.item()),'same pixel size in both lat and lon in xarray'
    
    #---- calculate the grid of test points and their weights
    dlat=dlon=pix_lat.item()
    xlat = []; xlon = []; area = []
    nlat = 180./dlat; nlon = 360./dlon
    if not np.mod(nlat,1.)==0.: raise AssertionError('nlat should be an integer')
    if not np.mod(nlon,1.)==0.: raise AssertionError('nlon should be an integer')
    nlat = int(nlat); nlon = int(nlon)
    for ilat in range(nlat):
        xlat.append(90.0-0.5*dlat-(ilat*dlat))
    for ilon in range(nlon):
        val = 0.5*dlon+(ilon*dlon)
        if min(data.coords[lonname].values) < 0.: # if xarray goes from (-180,180)
            if val>360.:
                xlon.append(val-360.)
            else:
                xlon.append(val)
        else:
            xlon.append(val)
    for ilat in range(nlat):
        area.append(2.*np.pi*(sind(xlat[ilat]+0.5*dlat)-             sind(xlat[ilat]-0.5*dlat))/float(nlon))
    
    # now fill the areas
    areaarray = []
    if data.shape[0] == 0.5*data.shape[1]:
        rowvar = latname
        colvar = lonname
        latdim = 'row'
    elif data.shape[1] == 0.5*data.shape[0]:
        rowvar = lonname
        colvar = latname
        latdim = 'col'
    else:
        raise ValueError('dimensions should be data.shape[0] == 0.5*data.shape[1]')
    areaarray = np.zeros(data.shape)        
    if latdim == 'row':
        for irow in range(len(data.coords[rowvar])):
            ifind = int((90.0-0.5*dlat-data.coords[rowvar][irow].item())/dlat)
            areaarray[ifind,:]=area[ifind]
    elif latdim == 'col':
        for icol in range(len(data.coords[colvar])):
            ifind = int((90.0-0.5*dlat-data.coords[colvar][icol].item())/dlat)
            areaarray[:,ifind]=area[ifind]
        
    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    area  = xr.DataArray(areaarray,name='area',coords=data.drop(drops).coords)
    return area


def MeanDataArray(data,area=None,latname = 'latitude', lonname = 'longitude'):
    """
    weighted average for xray data geographically averaged

    Parameters
    ----------
    data : DataArray
        the xray object to average over
        
    area: DataArray containing area weights. Obtained from AreaDataArray. 
          If None, calculate again.

    Returns
    -------
    
    globalav: global average

    area : DataArray containing area weights
    
    percentglobal: percentage of global area covered by this basis set
    
    """

    if not isinstance(data, xr.DataArray): raise ValueError("date must be an xray DataArray")
    
    pix_lat = np.unique(np.ediff1d(np.sort(data.coords[latname].values)))
    pix_lon = np.unique(np.ediff1d(np.sort(data.coords[lonname].values)))
    if not len(pix_lat)==len(pix_lon)==1: raise AssertionError('only one pixel size allowed in xarray')
    if not pix_lat.item()==pix_lon.item(): raise AssertionError('same pixel size in both lat and lon in xarray')
    # take weights
    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    try:
        totarea = np.sum(area.values)
    except:# if area does not exist, evaluate it
        area = AreaDataArray(data.drop(drops),latname,lonname)
        totarea = np.sum(area.values)
    percentglobal = np.round(totarea/(4.*np.pi)*100.,3)
    weighted = area*data
    # find the precision of data to round the average to
    max_precision = 0
    for val in data.drop(drops).values.flatten():
        _ , precision = precision_and_scale(val)
        if precision > max_precision: max_precision = precision
    average = np.round(np.sum(weighted.values)/totarea,decimals=max_precision)
    return average,area,percentglobal
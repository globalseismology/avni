#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int

import os
import numpy as np
import xarray as xr
from scipy.spatial import cKDTree
import pickle
import warnings
import pdb

####################### IMPORT REM3D LIBRARIES  #######################################
from .trigd import sind
from ..mapping import spher2cart
from .. import constants
from .common import precision_and_scale,convert2nparray
from ..tools import decimals
#######################################################################################

def tree3D(treefile,latitude=None,longitude=None,radius_in_km=None):
    #Build the tree if none is provided
    if os.path.isfile(treefile):
        print('... Reading KDtree file '+treefile)
        tree = pickle.load(open(treefile,'rb'))
    else:
        print('... KDtree file '+treefile+' not found for interpolations. Building it')
        if np.any(latitude == None) or np.any(longitude == None) or np.any(radius_in_km == None):
            raise IOError('latitude, longitude and radius_in_km are required')
        latitude = convert2nparray(latitude)
        longitude = convert2nparray(longitude)
        radius_in_km = convert2nparray(radius_in_km)
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

def ncfile2tree3D(ncfile,treefile,lonlatdepth = None,stride=None, radius_in_km = None):
    """
    Read or write a pickle interpolant with KDTree

    Input Parameters:
    ----------------

    ncfile: netcdf format file e.g. topography

    stride: is the downsampling before interpolation.

    radius_in_km: radius in kilometer when a 2D surface is valid. Ignores the
                    3rd field in lonlatdepth.Typically 6371km for Earth

    lonlatdepth: variable name of the longitude, latitude, depth (in km) arrays
                 default: ['longitude','latitude','depth']
    """
    #defaults
    if lonlatdepth is None: lonlatdepth = ['longitude','latitude','depth']

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
        rad = constants.R.to('km').magnitude - dep
    else:
        rad = xr.IndexVariable('rad',[radius_in_km])
    f.close() #close netcdf file

    # get the tree
    gridlat, gridrad, gridlon = np.meshgrid(lat.data,rad.data,lon.data)
    tree = tree3D(treefile,gridlat,gridlon,gridrad)
    return tree

def checkDataArray(data,latname = 'latitude', lonname = 'longitude'):
    """
    checks whether the data input is a DataArray and the coordinates are compatible

    Parameters
    ----------
    data : Dataset or DataArray
        the xray object to average over

    Returns
    -------

    pix: size of the uniform pixel

    """
    if not isinstance(data, xr.DataArray): raise ValueError("date must be an xray DataArray")

    pix_lat = np.unique(np.ediff1d(np.sort(data.coords[latname].values)))
    pix_lon = np.unique(np.ediff1d(np.sort(data.coords[lonname].values)))
    # restrict comparison to the minimum number of decimal points
    # dec_lat = decimals(pix_lat)
    # dec_lon = decimals(pix_lon)
    # pix_lat = np.unique(np.round(pix_lat, min(dec_lat)))
    # pix_lon = np.unique(np.round(pix_lon, min(dec_lon)))

    # checks
    ierror = []
    if len(pix_lat)==len(pix_lon)==1:
        if not pix_lat.item()==pix_lon.item(): warnings.warn('same pixel size in both lat and lon in xarray')
        # check number of poxels
        pix_width = pix_lat.item()
        nlat = 180./pix_width; nlon = 360./pix_width
        if np.mod(nlat,1.) is not 0. or np.mod(nlon,1.) is not 0.:
            ierror.append(1)
            warnings.warn ('pixel width should be ideally be a factor of 180')
        nlat = int(nlat); nlon = int(nlon)
        if nlat*nlon != data.size:
            ierror.append(2)
            warnings.warn('number of pixels expected for '+str(pix_width)+'X'+str(pix_width)+' is '+str(nlat*nlon)+',  not '+str(data.size)+' as specified in the data array.')
    else:
        warnings.warn('multiple pixel sizes have been found for xarray'+str(pix_lat)+str(pix_lon)+'. Choosing the first value '+str(pix_lat[0]))
        ierror.append(3)
        pix_width = pix_lat[0]
    return ierror,pix_width


def AreaDataArray(data,latname = 'latitude', lonname = 'longitude',pix_width=None):
    """
    weighted average for xray data geographically averaged

    Parameters
    ----------
    data : Dataset or DataArray
        the xray object to average over

    pix_width: width of pixels if not the default derived from data

    Returns
    -------

    area: a DataArray object with area of each pixel

    """

    # check if it is a compatible dataarray
    ierror,pix = checkDataArray(data,latname, lonname)

    if pix_width is not None:
        if pix_width.shape != data.shape: raise AssertionError('pix_width.shape != data.shape')
        uniq_pix = np.unique(pix_width)
        # if the pix_width array has only one value and that
        # is consistent with the one derived from data
        if len(uniq_pix) is 1 and uniq_pix[0] is pix: pix_width = None

    # now fill the areas
    area = {}
    areaarray = np.zeros(data.shape)

    # find the index of the latitude
    lat_index = data.dims.index(latname)
    if lat_index is not 0: # transpose to (lat,lon) for caculations
        data = data.T
        if pix_width is not None: pix_width = pix_width.T

    # now calculate area
    for irow in range(len(data.coords[latname])):
        xlat = data.coords[latname][irow].item()
        # if no px width array is provided
        if pix_width is None:
            dlat = dlon = pix
            ifind = int((90.0-0.5*dlat-xlat)/dlat)
            if ifind not in area.keys():
                nlon = int(180./pix)
                area[ifind] = 2.*np.pi*(sind(xlat+0.5*dlat)-             sind(xlat-0.5*dlat))/float(nlon)
            areaarray[ifind,:]=area[ifind]
        else:
            for icol in range(len(pix_width[lonname])):
                nlon = int(180./pix_width[irow,icol])
                dlat = dlon = pix_width[irow,icol].item()
                areaarray[irow,icol] = 2.*np.pi*(sind(xlat+0.5*dlat)-             sind(xlat-0.5*dlat))/float(nlon)

    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    area  = xr.DataArray(areaarray,name='area',coords=data.drop(drops).coords)

    # transpose the area array if needed
    if lat_index is not 0: area = area.T
    return area


def MeanDataArray(data,area=None,latname = 'latitude', lonname = 'longitude', pix_width=None):
    """
    weighted average for xray data geographically averaged

    Parameters
    ----------
    data : DataArray
        the xray object to average over

    area: DataArray containing area weights. Obtained from AreaDataArray.
          If None, calculate again.

    pix_width: width of pixels if not the default derived from data

    Returns
    -------

    globalav: global average

    area : DataArray containing area weights

    percentglobal: percentage of global area covered by this basis set

    """
    if pix_width is not None:
        if pix_width.shape != data.shape: raise AssertionError('pix_width.shape != data.shape')

    # take weights
    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    if area is None: area = AreaDataArray(data.drop(drops),latname,lonname,pix_width)
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

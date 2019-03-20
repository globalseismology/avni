#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import numpy as np
import xarray as xr

####################### IMPORT REM3D LIBRARIES  #######################################
from .trigd import sind
#######################################################################################

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
    assert(np.mod(nlat,1.)==0.),'nlat should be an integer'
    assert(np.mod(nlon,1.)==0.),'nlon should be an integer'
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
    totarea=0.; areaarray = []
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
    for irow in range(len(data.coords[rowvar])):
        arearow = []
        for icol in range(len(data.coords[colvar])):
            if latdim == 'row':
                arearow.append(area[irow])
                totarea=totarea+area[irow]
            elif latdim == 'col':
                arearow.append(area[icol])
                totarea=totarea+area[icol]
        areaarray.append(arearow)
        
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
    assert(len(pix_lat)==len(pix_lon)==1),'only one pixel size allowed in xarray'
    assert(pix_lat.item()==pix_lon.item()),'same pixel size in both lat and lon in xarray'
    # take weights
    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    area  = AreaDataArray(data.drop(drops),latname,lonname)    
    totarea = np.sum(area.values)
    percentglobal = round(totarea/(4.*np.pi)*100.,3)
    weighted = area*data
    average = np.mean(weighted.values)    
    return average,area,percentglobal
#!/usr/bin/env python

#####################  IMPORT STANDARD MODULES   #########################

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
import math
from scipy import sparse
import typing as tp
#from mpl_toolkits.basemap import shiftgrid

####################### IMPORT AVNI LIBRARIES  ###########################

from .trigd import sind
from ..mapping import spher2cart
from .. import constants
from .common import precision_and_scale,convert2nparray
from ..tools import get_filedir
from ..data import common

##########################################################################

def xarray_to_epix(data: tp.Union[xr.DataArray,xr.Dataset],
                   latname: str = 'latitude',
                   lonname: str = 'longitude') -> np.ndarray:
    """Convert multi-dimensional pixel grid in xarray to extended pixel (.epix) format

    Parameters
    ----------
    data : tp.Union[xr.DataArray,xr.Dataset]
        Multi-dimensional pixel grid in xarray formats
    latname : str, optional
        Name to use for latitude column in named numpy array, by default 'latitude'
    lonname : str, optional
        Name to use for longitude column in named numpy array, by default 'longitude'

    Returns
    -------
    np.ndarray
        Array containing (`latitude`, `longitude`, `pixel_size`, `value`)

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # check if it is a compatible xarray
    _,pix,_ = checkxarray(data,latname, lonname)

    if data.dims[0] == latname:
        values = data.T.data.ravel()
        nlat = data.shape[0]; nlon = data.shape[1]
    else:
        values = data.data.ravel()
        nlat = data.shape[1]; nlon = data.shape[0]
    lons = np.repeat(data[lonname].data,nlat)
    lats = np.tile(data[latname].data,nlon)
    pixsize = pix*np.ones_like(lons)
    epixarr = np.vstack((lats,lons,pixsize,values)).T
    dt = {'names':[latname, lonname,'pixel_size','value'], 'formats':[float, float,float,float]}
    epixarr = np.zeros(len(lats),dtype=dt)
    epixarr[latname] = lats
    epixarr[lonname] = lons
    epixarr['pixel_size'] = pixsize
    epixarr['value'] = values
    return epixarr

def epix_to_xarray(epixarr: np.ndarray,
                   latname: str = 'latitude',
                   lonname: str = 'longitude') -> xr.DataArray:
    """Convert extended pixel (.epix) to multi-dimensional pixel grid in xarray format

    Parameters
    ----------
    epixarr : np.ndarray
        Array containing (`latitude`, `longitude`, `pixel_size`, `value`)
    latname : str, optional
        Name to use for latitude column in named numpy array, by default 'latitude'
    lonname : str, optional
        Name to use for longitude column in named numpy array, by default 'longitude'

    Returns
    -------
    xr.DataArray
        Multi-dimensional pixel grid in xarray formats

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    lonspace = np.setdiff1d(np.unique(np.ediff1d(np.sort(epixarr[lonname]))),[0.])
    latspace = np.setdiff1d(np.unique(np.ediff1d(np.sort(epixarr[latname]))),[0.])
    if not(len(lonspace) == len(lonspace) == 1): raise AssertionError('not(len(lonspace) == len(lonspace) == 1)')
    spacing = lonspace[0]
    latitude = np.arange(-90.+spacing/2.,90.,spacing)
    longitude = np.arange(0.+spacing/2.,360.,spacing)
    outarr = xr.DataArray(np.zeros((len(latitude),len(longitude))),
                    dims=[latname, lonname],
                    coords={latname:latitude,lonname:longitude})
    for ii, lat in enumerate(epixarr[latname]):
        lon = epixarr[lonname][ii]
        val = epixarr['value'][ii]
        if spacing != epixarr['pixel_size'][ii]: raise ValueError('spacing != epixarr[pixel_size]')
        outarr.loc[dict(latitude=lat,longitude=lon)] = val
    return outarr

def tree3D(treefile: str,
           latitude: tp.Union[None,list,tuple,np.ndarray] = None,
           longitude: tp.Union[None,list,tuple,np.ndarray] = None,
           radius_in_km: tp.Union[None,list,tuple,np.ndarray] = None):
    """Build a KD-tree at specific locations

    Parameters
    ----------
    treefile : str
        Name of the file where tree is (or will be) stored
    latitude : tp.Union[list,tuple,np.ndarray], optional
        Latitudes of locations queried, by default None
    longitude : tp.Union[list,tuple,np.ndarray], optional
        Longitudes of locations queried, by default None
    radius_in_km : tp.Union[list,tuple,np.ndarray], optional
        Radii of locations queried, by default None

    Returns
    -------
    tree
        A :py:func:`scipy.spatial.cKDTree`

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    #Build the tree if none is provided
    build_tree = False
    if os.path.isfile(treefile):
        print('... Reading KDtree file '+treefile)
        try:
            tree = pickle.load(open(treefile,'rb'))
        except:
            build_tree = True
    else:
        build_tree = True

    if build_tree:
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

def querytree3D(tree,
                latitude: tp.Union[list,tuple,np.ndarray],
                longitude: tp.Union[list,tuple,np.ndarray],
                radius_in_km: tp.Union[list,tuple,np.ndarray],
                values: tp.Union[None,sparse.csc_matrix,sparse.csr_matrix,list,tuple,np.ndarray] = None,
                nearest: int = 1):
    """Query a KD-tree for values at specific locations

    Parameters
    ----------
    tree
        A :py:func:`scipy.spatial.cKDTree`
    latitude : tp.Union[list,tuple,np.ndarray], optional
        Latitudes of locations queried, by default None
    longitude : tp.Union[list,tuple,np.ndarray], optional
        Longitudes of locations queried, by default None
    radius_in_km : tp.Union[list,tuple,np.ndarray], optional
        Radii of locations queried, by default None
    values : tp.Union[None,sparse.csc_matrix,sparse.csr_matrix,list,tuple,np.ndarray], optional
        Values at each KD-tree point, by default None
    nearest : int, optional
        Number of nearest values in the KD-tree to interpolated from, by default
        1 so nearest

    Returns
    -------
    inds or interp, inds
        indices of the nearest points in te KD-tree and the interporlated value
        (if values is not None)

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # checks
    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)
    radius_in_km = convert2nparray(radius_in_km)
    if not(len(latitude)==len(latitude)==len(radius_in_km)): raise AssertionError('latitude, longitude and radius need to be of same size')
    evalpoints = np.column_stack((radius_in_km.flatten(),latitude.flatten(), longitude.flatten()))
    coordstack = spher2cart(evalpoints)
    dist,inds = tree.query(coordstack,k=nearest)
    if np.any(values==None):
        return inds
    else:
        # convert to csc_matrix
        if not isinstance(values,sparse.csc_matrix):
            if isinstance(values,sparse.csr_matrix):
                values = values.tocsc()
            else:
                # do not allow strings as values
                values = convert2nparray(values,int2float = True,allowstrings=False)
                if values.ndim != 1: raise ValueError('only 1 dimenional values allowed')
                values = sparse.csc_matrix(values).transpose().tocsc()

        # find values
        if nearest == 1:
            interp = values[inds]
        else:
            weights = 1.0 / dist**2
            rows = ((np.arange(inds.shape[0])*np.ones_like(inds).T).T).ravel()
            cols = inds.ravel()
            weighted = values[cols].reshape(inds.shape,order='C').multiply(weights)
            weighted_sum = np.asarray(weighted.sum(axis=1))
            weights_sum = weights.sum(axis=1).reshape(weighted_sum.shape)
            val_temp = np.divide(weighted_sum,weights_sum,out=np.zeros_like(weights_sum), where=weights_sum!=0)
            interp = sparse.csc_matrix(val_temp)
        return interp,inds

def get_stride(resolution: str) -> int:
    """Get the stride to use for various resolution of plotting. This dictates
    downsampling before interpolation.

    Parameters
    ----------
    resolution : str
        Resolution of boundary database to use in Basemap.
        Can be c (crude), l (low), i (intermediate), h (high), f (full)

    Returns
    -------
    int
        stride to use for resolution

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    if resolution.lower() == 'c' or resolution.lower() == 'crude':
        stride = 20
    elif resolution.lower() == 'l' or resolution.lower() == 'low':
        stride = 10
    elif resolution.lower() == 'i' or resolution.lower() == 'intermediate':
        stride = 5
    elif resolution.lower() == 'h' or resolution.lower() == 'high':
        stride = 2
    elif resolution.lower() == 'f' or resolution.lower() == 'full':
        stride = 1
    else:
        raise KeyError('resolution can only be low, medium, high')
    return stride

def ncfile2tree3D(ncfile: str,treefile: str,
                  lonlatdepth: list = ['longitude','latitude','depth'],
                  resolution: str = 'h',
                  radius_in_km: tp.Union[None,list,tuple,np.ndarray] = None):
    """Read or write a pickle interpolant with KD-tree

    Parameters
    ----------
    ncfile : str
        Name of the topography file in NETCDF4 format
    treefile : str
        Name of the file where tree is (or will be) stored
    lonlatdepth : list, optional
        A list of variable names of the longitude, latitude, depth (in km) arrays, by default ['longitude','latitude','depth']
    resolution : str, optional
        Dictates downsampling before interpolation, by default 'h'
    radius_in_km : tp.Union[None,list,tuple,np.ndarray], optional
        Radius in kilometer when a 2D surface. Ignores the 3rd field in `lonlatdepth`, by default None

    Returns
    -------
    tree
        A :py:func:`scipy.spatial.cKDTree`

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    #read topography file
    stride = get_stride(resolution)
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

def readtopography(model: tp.Union[None,str] = None,
                   resolution: str = 'h',
                   field: str = 'z',
                   latitude: str = 'lat',
                   longitude: str = 'lon',
                   latitude_limits: list = [-90,90],
                   longitude_limits: list = [-180,180],
                   dbs_path: tp.Union[None,str] = None):
    """Read standard topography file in NETCDF4 format specified in :py:func:`constants`

    Parameters
    ----------
    model : str, optional
        Name of the topography file in NETCDF4 format, by default :py:func:`constants.topography`
    resolution : str, optional
        Dictates downsampling before interpolation, by default 'h'
    field : str, optional
        Field name in the NETCDF4 file to use, by default 'z'
    latitude : str, optional
        Name to use for latitude column in named numpy array, by default 'lat'
    longitude : str, optional
        Name to use for longitude column in named numpy array, by default 'lon'
    latitude_limits : list, optional
        Limit for restricting the domain for reading topography, by default [-90,90]
    longitude_limits : list, optional
        Limit for restricting the domain for reading topography, by default [-180,180]
    dbs_path : tp.Union[None,str], optional
        Database path to the folder contianing topography file, by default None so takes the value in :py:func:`constants.topofolder`

    Returns
    -------
    xr.Dataset
        Multi-dimensional pixel grid in xarray formats

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # Get the directory location where CPT files are kep
    if dbs_path is None: dbs_path = get_filedir(subdirectory=constants.topofolder)
    if model is None: model = constants.topography

    # Download file if possible
    ncfile = os.path.join(dbs_path,model)
    if not os.path.isfile(ncfile):
        success = False
        _,success = common.update_file(model,subdirectory=constants.topofolder)
        if not success: ValueError("Could not find file "+model)

    f = xr.open_dataset(ncfile)

    # Get the stride based on requested resolution
    stride = get_stride(resolution)
    model = f[field][::stride,::stride]

    # subselect region within -not implemented yet
    #shift it by the grid
#     valout, lonout = shiftgrid(longitude_limits[0],model.data,model['lon'].data,start=True)
#     lonout[np.where(lonout>180)] =lonout[lonout>180]-360.
#     model[longitude].data = lonout
#     model.data  = valout
    model  = model.sortby('lon') # required for transform_scalar
    f.close() #close netcdf file
    return model

def checkxarray(data: tp.Union[xr.DataArray,xr.Dataset],
                latname: str = 'latitude',
                lonname: str = 'longitude',
                Dataset = True) -> tp.Tuple[list,float,tuple]:
    """Checks whether the data input is a DataArray and the coordinates are compatible

    Parameters
    ----------
    data : tp.Union[xr.DataArray,xr.Dataset]
        Multi-dimensional pixel grid in xarray formats
    latname : str, optional
        Name to use for latitude column in named numpy array, by default 'latitude'
    lonname : str, optional
        Name to use for longitude column in named numpy array, by default 'longitude'
    Dataset : bool, optional
        Allow Dataset or not, by default True

    Returns
    -------
    tp.Tuple[list,float,tuple]
        First element is a list of error warnings while performing checks

        Second element is size of the uniform pixel.

        Third element is a tuple containing the shape of the grid

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    if not isinstance(data, (xr.DataArray,xr.Dataset)): raise ValueError("date must be an xarray DataArray or Dataset")
    if not Dataset and isinstance(data, xr.Dataset): raise ValueError("date must be an xarray Dataset if True")

    pix_lat = np.unique(np.ediff1d(np.sort(data.coords[latname].values)))
    pix_lon = np.unique(np.ediff1d(np.sort(data.coords[lonname].values)))
    if isinstance(data, xr.DataArray) : shape = data.shape
    if isinstance(data, xr.Dataset) : shape = tuple(data.dims[d] for d in [latname, lonname])

    # checks
    ierror = []
    if len(pix_lat)==len(pix_lon)==1:
        if not pix_lat.item()==pix_lon.item(): warnings.warn('same pixel size in both lat and lon in xarray')
        # check number of poxels
        pix_width = pix_lat.item()
        if not math.isclose(180 % pix_width,0):
            ierror.append(1)
            warnings.warn ('pixel width should be ideally be a factor of 180')
        nlat = int(180./pix_width); nlon = int(360./pix_width)
        if nlat*nlon != shape[0]*shape[1]:
            ierror.append(2)
            warnings.warn('number of pixels expected for '+str(pix_width)+'X'+str(pix_width)+' is '+str(nlat*nlon)+',  not '+str(data.size)+' as specified in the data array.')
    else:
        warnings.warn('multiple pixel sizes have been found for xarray'+str(pix_lat)+str(pix_lon)+'. Choosing the first value '+str(pix_lat[0]))
        ierror.append(3)
        pix_width = pix_lat[0]
    return ierror,pix_width,shape


def areaxarray(data: tp.Union[xr.DataArray,xr.Dataset],
               latname: str = 'latitude',
               lonname: str = 'longitude',
               pix_width: tp.Union[None,np.ndarray] = None) -> xr.DataArray:
    """Calculate area for multi-dimensional pixel grid in xarray formats

    Parameters
    ----------
    data : tp.Union[xr.DataArray,xr.Dataset]
        Multi-dimensional pixel grid in xarray formats
    latname : str, optional
        Name to use for latitude column, by default 'latitude'
    lonname : str, optional
        Name to use for longitude column, by default 'longitude'
    pix_width : tp.Union[None,np.ndarray], optional
        Width of pixels if not the default derived from data, by default None so
        is derived

    Returns
    -------
    xr.DataArray
        A DataArray object with area of each pixel

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # check if it is a compatible dataarray
    ierror,pix,shape = checkxarray(data,latname, lonname)

    if pix_width is not None:
        if pix_width.shape != shape: raise AssertionError('pix_width.shape != shape')
        uniq_pix = np.unique(pix_width)
        # if the pix_width array has only one value and that
        # is consistent with the one derived from data
        if len(uniq_pix) == 1 and uniq_pix[0] == pix: pix_width = None

    # now fill the areas
    area = {}
    areaarray = np.zeros(shape)

    # find the index of the latitude
    lat_index = np.argwhere(np.array(data.dims)==latname)[0].item()
    if lat_index != 0: # transpose to (lat,lon) for caculations
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
                nlon = int(360./pix)
                area[ifind] = 2.*np.pi*(sind(xlat+0.5*dlat)-sind(xlat-0.5*dlat))/float(nlon)
            areaarray[ifind,:]=area[ifind]
        else:
            for icol in range(len(pix_width[lonname])):
                nlon = int(360./pix_width[irow,icol])
                dlat = dlon = pix_width[irow,icol].item()
                areaarray[irow,icol] = 2.*np.pi*(sind(xlat+0.5*dlat)-sind(xlat-0.5*dlat))/float(nlon)

    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    area  = xr.DataArray(areaarray,name='area',coords=data.drop(drops).coords)

    # transpose the area array if needed
    if lat_index != 0: area = area.T
    return area

def meanxarray(data: tp.Union[xr.DataArray,xr.Dataset],
               area: tp.Union[None,xr.DataArray] = None,
               latname: str = 'latitude',
               lonname: str = 'longitude',
               pix_width: tp.Union[None,np.ndarray] = None) -> tp.Tuple[float,xr.DataArray,float]:
    """Calculate geographically weighted average of a multi-dimensional pixel grid in xarray formats

    Parameters
    ----------
    data : tp.Union[xr.DataArray,xr.Dataset]
        Multi-dimensional pixel grid in xarray formats to average over
    area : tp.Union[None,xr.DataArray], optional
        Area of each pixel, by default None so calculated on the fly
    latname : str, optional
        Name to use for latitude column, by default 'latitude'
    lonname : str, optional
        Name to use for longitude column, by default 'longitude'
    pix_width : tp.Union[None,np.ndarray], optional
        Width of pixels if not the default derived from data, by default None so
        is derived

    Returns
    -------
    tp.Tuple[float,xr.DataArray,float]
        First element is the global average
        Second element is a DataArray containing area weights for each pixel
        Third element is the percentage of global area covered by this basis set

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # check if it is a compatible dataarray
    ierror,pix,shape = checkxarray(data,latname, lonname)

    if pix_width is not None:
        if pix_width.shape != shape: raise AssertionError('pix_width.shape != shape')

    # take weights
    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    if area is None: area = areaxarray(data.drop(drops),latname,lonname,pix_width)
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
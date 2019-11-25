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
import math
from scipy import sparse
from mpl_toolkits.basemap import shiftgrid

####################### IMPORT REM3D LIBRARIES  #######################################
from .trigd import sind
from ..mapping import spher2cart,intersection,midpoint
from .. import constants
from .common import precision_and_scale,convert2nparray
from ..tools import decimals,get_filedir
#######################################################################################

def xarray_to_epix(data,latname = 'latitude', lonname = 'longitude'):
    # check if it is a compatible dataarray
    ierror,pix,shape = checkxarray(data,latname, lonname)

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
    dt = {'names':[latname, lonname,'pixel_size','value'], 'formats':[np.float, np.float,np.float,np.float]}
    epixarr = np.zeros(len(lats),dtype=dt)
    epixarr[latname] = lats
    epixarr[lonname] = lons
    epixarr['pixel_size'] = pixsize
    epixarr['value'] = values
    return epixarr

def epix_to_xarray(epixarr,latname = 'latitude', lonname = 'longitude'):
    lonspace = np.setdiff1d(np.unique(np.ediff1d(np.sort(epixarr[lonname]))),[0.])
    latspace = np.setdiff1d(np.unique(np.ediff1d(np.sort(epixarr[latname]))),[0.])
    if not(len(lonspace) == len(lonspace) == 1): raise AssertionError('not(len(lonspace) == len(lonspace) == 1)')
    spacing = lonspace[0]
    latitude = np.arange(-90.+spacing/2.,90.,1.)
    longitude = np.arange(0.+spacing/2.,360.,1.)
    outarr = xr.DataArray(np.zeros((len(latitude),len(longitude))),
                    dims=[latname, lonname],
                    coords={latname:latitude,lonname:longitude})
    for ii, lat in enumerate(epixarr[latname]):
        lon = epixarr[lonname][ii]
        val = epixarr['value'][ii]
        if spacing != epixarr['pixel_size'][ii]: raise ValueError('spacing != epixarr[pixel_size]')
        outarr.loc[dict(latitude=lat,longitude=lon)] = val
    return outarr

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

def querytree3D(tree,latitude,longitude,radius_in_km,values=None,nearest=1):
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

def get_stride(resolution):
    """
    resolution: of boundary database to use like in Basemap.
                Can be c (crude), l (low), i (intermediate), h (high), f (full)
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

def ncfile2tree3D(ncfile,treefile,lonlatdepth = None,resolution='h', radius_in_km = None):
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

def readtopography(model=constants.topography,resolution='h',field = 'z',latitude='lat',longitude='lon',latitude_limits=[-90,90],longitude_limits=[-180,180], dbs_path=get_filedir()):
    """
    resolution: stride is 1 for high, 10 for low
    """
    stride = get_stride(resolution)
    if dbs_path == None:
        ncfile = model
    else:
        ncfile = dbs_path+'/'+model
    if not os.path.isfile(ncfile): data.update_file(model)
    #read values
    if os.path.isfile(ncfile):
        f = xr.open_dataset(ncfile)
    else:
        raise ValueError("Error: Could not find file "+ncfile)

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

def checkxarray(data,latname = 'latitude', lonname = 'longitude', Dataset = True):
    """
    checks whether the data input is a DataArray and the coordinates are compatible

    Parameters
    ----------
    data : Dataset or DataArray
        the xray object to average over

    Dataset : allow Dataset or not

    Returns
    -------

    pix: size of the uniform pixel

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


def areaxarray(data,latname = 'latitude', lonname = 'longitude',pix_width=None):
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
    ierror,pix,shape = checkxarray(data,latname, lonname)

    if pix_width is not None:
        if pix_width.shape != shape: raise AssertionError('pix_width.shape != shape')
        uniq_pix = np.unique(pix_width)
        # if the pix_width array has only one value and that
        # is consistent with the one derived from data
        if len(uniq_pix) is 1 and uniq_pix[0] is pix: pix_width = None

    # now fill the areas
    area = {}
    areaarray = np.zeros(shape)

    # find the index of the latitude
    lat_index = np.argwhere(np.array(data.dims)==latname)[0].item()
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
                nlon = int(360./pix)
                area[ifind] = 2.*np.pi*(sind(xlat+0.5*dlat)-             sind(xlat-0.5*dlat))/float(nlon)
            areaarray[ifind,:]=area[ifind]
        else:
            for icol in range(len(pix_width[lonname])):
                nlon = int(360./pix_width[irow,icol])
                dlat = dlon = pix_width[irow,icol].item()
                areaarray[irow,icol] = 2.*np.pi*(sind(xlat+0.5*dlat)-             sind(xlat-0.5*dlat))/float(nlon)

    # drop the variables for weights
    drops = [var for var in data.coords.keys() if var not in [latname,lonname]]
    area  = xr.DataArray(areaarray,name='area',coords=data.drop(drops).coords)

    # transpose the area array if needed
    if lat_index is not 0: area = area.T
    return area

def meanxarray(data,area=None,latname = 'latitude', lonname = 'longitude', pix_width=None):
    """
    weighted average for xray data geographically averaged

    Parameters
    ----------
    data : DataArray
        the xray object to average over

    area: DataArray containing area weights. Obtained from areaxarray.
          If None, calculate again.

    pix_width: width of pixels if not the default derived from data

    Returns
    -------

    globalav: global average

    area : DataArray containing area weights

    percentglobal: percentage of global area covered by this basis set

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

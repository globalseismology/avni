#!/usr/bin/env python

#####################  IMPORT STANDARD MODULES   #########################

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import numpy as np
from collections import Counter
from scipy import sparse
from timeit import default_timer as timer
from numba import jit,int64
import typing as tp

####################### IMPORT AVNI LIBRARIES  ###########################

from avni.f2py import vbspl,dbsplrem,ylm,shold
from .trigd import sind,cosd,acosd
from .common import convert2nparray,makegrid

##########################################################################

def eval_vbspl(radius: tp.Union[list,tuple,np.ndarray],
               knots: tp.Union[str,list,tuple,float,np.int64,np.ndarray,bool]):
    """Evaluate the cubic spline knot with second derivative as 0 at end points.
       In contrast to :py:func:`eval_splrem`, this function distributes the
       spline knots unevenly at the specified depths or radii.

    Parameters
    ----------
    radius : tp.Union[list,tuple,np.ndarray]
        A radius value or array of radii (or alternatively depths) queried in km
    knots : tp.Union[str,list,tuple,float,np.int64,np.ndarray,bool]
        A numpy array or list of radii (or depths) of spline knots in km

    Returns
    -------
    vercof, dvercof: float or np.ndarray
        value of the polynomial coefficients at each depth and derivative.
        Both arrays have size (Nradius, Nsplines).

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    if isinstance(knots, (list,tuple,np.ndarray)):
        knots = np.asarray(knots)
        knots = np.sort(knots)
    else:
        raise TypeError('knots must be list or tuple, not %s' % type(knots))

    # convert to numpy arrays
    radius = convert2nparray(radius)

    # find repeated values
    repeats = [item for item, count in Counter(knots).items() if count > 1]
    repeats_gt_2= [item for item, count in Counter(knots).items() if count > 2]
    if len(repeats_gt_2) != 0: raise ValueError('Cannot have more than 2 repetitions in knots')

    # if there are repeated knots, splits it
    if len(repeats) > 0:
        split_at = []
        for index,value in enumerate(repeats):
            split_at.append(np.where(knots==value)[0][1])
        knots_list = np.split(knots, split_at)
        for knots in knots_list:
            if len(knots) < 4:
                raise ValueError('Atleast 4 knots need to be defined at or between '+str(min(knots))+' and '+str(max(knots))+' km')

        jj = 0
        for query in radius:
            jj = jj + 1
            for index,value in enumerate(knots_list):
                # create the arrays as Fortran-contiguous
                splpts = np.array(value.tolist(), order = 'F')
                #Undefined if query does not lie within the query extents of knot points
                if query < min(value) or query > max(value):
                    temp1 = temp2 = np.zeros_like(splpts)
                else:
                    (temp1, temp2) = vbspl(query,len(splpts),splpts)
                if index == 0:
                    vercof_temp = temp1; dvercof_temp = temp2
                else:
                    vercof_temp = np.concatenate((vercof_temp,temp1))
                    dvercof_temp = np.concatenate((dvercof_temp,temp1))
            if jj == 1:
                vercof = vercof_temp; dvercof = dvercof_temp
            else:
                vercof = np.vstack([vercof,vercof_temp])
                dvercof = np.vstack([dvercof,dvercof_temp])

    # when there are no repeated knots
    else:
        if len(knots) < 4:
            raise ValueError('Atleast 4 knots need to be defined at or between '+str(min(knots))+' and '+str(max(knots)))
        # create the arrays as Fortran-contiguous
        splpts = np.array(knots.tolist(), order = 'F')
        # initialize
        vercof = np.ones((len(radius),len(knots)))
        dvercof = np.zeros((len(radius),len(knots)))
        for idep,query in enumerate(radius):
            #Undefined if query does not lie within the query extents of knot points
            if query < min(knots) or query > max(knots):
                vercof_temp = dvercof_temp = np.zeros_like(splpts)
            else:
                (vercof_temp, dvercof_temp) = vbspl(query,len(splpts),splpts)
            vercof[idep]=vercof_temp
            dvercof[idep]=dvercof_temp

    return vercof, dvercof


def eval_splrem(radius: tp.Union[list,tuple,np.ndarray],
                radius_range: tp.Union[list,tuple,np.ndarray],
                nsplines: int):
    """Evaluate the cubic spline knot with second derivative as 0 at end points.
       In contrast to :py:func:`eval_vbspl`, this function distributes the
       spline knots evenly within the radius range.


    Parameters
    ----------
    radius : tp.Union[list,tuple,np.ndarray]
        A radius value or array of radii (or alternatively depths) queried in km
    radius_range : tp.Union[list,tuple,np.ndarray]
        Limits of the radius (or depths) limits of the region
    nsplines : int
        number of splines within the range

    Returns
    -------
    vercof, dvercof: float or np.ndarray
        value of the polynomial coefficients at each depth and derivative.
        Both arrays have size (Nradius, Nsplines).

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # convert to numpy arrays
    radius = convert2nparray(radius)

    if len(radius_range) != 2 or not isinstance(radius_range, (list,tuple,np.ndarray)):
        raise TypeError('radius_range must be list , not %s' % type(radius_range))

    for irad,radval in enumerate(radius):
        #Undefined if depth does not lie within the depth extents of knot points
        if radval < min(radius_range) or radval > max(radius_range):
            temp1 = temp2 = np.zeros(nsplines)
        else:
            (temp1, temp2) = dbsplrem(radval,radius_range[0], radius_range[1],nsplines)
        if irad == 0:
            vercof = temp1; dvercof = temp2
        else:
            vercof = np.vstack([vercof,temp1])
            dvercof = np.vstack([dvercof,temp2])

    return vercof, dvercof


def eval_polynomial(radius: tp.Union[list,tuple,np.ndarray],
                    radius_range: tp.Union[list,tuple,np.ndarray],
                    rnorm: float,
                    types: list = ['CONSTANT','LINEAR']):
    """Evaluate a set of polynomial functions within a range of radii.

    Parameters
    ----------
    radius : tp.Union[list,tuple,np.ndarray]
        A radius value or array of radii queried in km
    radius_range : tp.Union[list,tuple,np.ndarray]
        Limits of the radius limits of the region
    rnorm : float
        normalization for radius, usually the radius of the planet
    types : list, optional
        polynomial coefficients to be used for calculation.
        Options are : TOP,BOTTOM, CONSTANT, LINEAR, QUADRATIC, CUBIC, by default ['CONSTANT','LINEAR']

    Returns
    -------
    vercof, dvercof: float or np.ndarray
        value of the polynomial coefficients and derivative at each depth size (Nradius)

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # convert to numpy arrays
    radiusin = convert2nparray(radius)
    radius_temp = convert2nparray(radius_range)
    if radius_temp.ndim == 1:
        radius_range = convert2nparray([radius_temp.tolist()])
        if radius_range.shape[1] != 2: raise TypeError('Only two values allowed within radius_range')
    elif radius_temp.ndim == 2:
        radius_range = radius_temp
        if radius_range.shape[1] != 2: raise TypeError('Only two values allowed within radius_range')
    if not isinstance(radius_range, (list,tuple,np.ndarray)):
        raise TypeError('radius_range must be list , not %s' % type(radius_range))

    # keys in coefficients should be acceptable
    choices = ['TOP', 'BOTTOM', 'CONSTANT', 'LINEAR', 'QUADRATIC', 'CUBIC']
    if not np.all([key in choices for key in types]): raise AssertionError('Only polynomial bases can be used')
    npoly = len(radius_range)*len(types)
    # first find whether CONSTANT/LINEAR or TOP/BOTTOM
    for radii in radius_range:
        if not np.all(np.sort(radii)==radii): raise AssertionError('radius_range needs to be sorted')

    # see if either TOP/BOT or CONSTANT/LINEAR exists
    findtopbot = np.any([key in ['BOTTOM','TOP'] for key in types])
    findconstantlinear = np.any([key in ['CONSTANT','LINEAR'] for key in types])
    if findtopbot and findconstantlinear: raise ValueError('Cannot have both BOTTOM/TOP and CONSTANT/LINEAR as types in eval_polynomial')

    # loop through each depth and get the polynomial values and derivatives
    for irad,_ in enumerate(radiusin):
        temp = np.zeros(npoly)
        dtemp = np.zeros(npoly)
        for irange,_ in enumerate(radius_range):
            rbot = min(radius_range[irange])/rnorm
            rtop = max(radius_range[irange])/rnorm
            #Undefined if depth does not lie within the depth extents of knot points
            if radiusin[irad]/rnorm < rbot or radiusin[irad]/rnorm > rtop:
                # <= so that the boundary depth belongs to only one radial kernel
                for itype in range(len(types)):
                    ii = irange*len(types)+itype
                    temp[ii]=0.
                    dtemp[ii]=0.
            else:
                rn=radiusin[irad]/rnorm
                for itype,_ in enumerate(types):
                    ii = irange*len(types)+itype
                    if findconstantlinear:
                        if types[itype]=='CONSTANT':
                            temp[ii]=1.
                            dtemp[ii]=0.
                        elif types[itype]=='LINEAR':
                            temp[ii]=rn
                            dtemp[ii]=1.
                        elif types[itype]=='QUADRATIC':
                            temp[ii]=rn**2
                            dtemp[ii]=2.*rn
                        elif types[itype]=='CUBIC':
                            temp[ii]=rn**3
                            dtemp[ii]=3.*rn**2
                    elif findtopbot:
                        if types[itype]=='TOP':
                            temp[ii] = 1.-(rn-rtop)/(rbot-rtop)
                            dtemp[ii]= -1./(rbot-rtop)
                        elif types[itype]=='BOTTOM':
                            temp[ii] = (rn-rtop)/(rbot-rtop)
                            dtemp[ii]= 1./(rbot-rtop)
                        elif types[itype]=='QUADRATIC':
                            temp[ii] = rn**2-rtop**2-(rn-rtop)*(rbot+rtop)
                            dtemp[ii]=2.*rn - 1./(rbot+rtop)
                        elif types[itype]=='CUBIC':
                            temp[ii]=rn**3-rtop**3-(rn-rtop)*(rbot**3-rtop**3)/(rbot-rtop)
                            dtemp[ii]=3.*rn**2 - 1.*(rbot**3-rtop**3)/(rbot-rtop)
        if irad == 0:
            vercof = temp; dvercof = dtemp
        else:
            vercof = np.vstack([vercof,temp])
            dvercof = np.vstack([dvercof,dtemp])

    return vercof,dvercof

def eval_splcon(latitude: tp.Union[list,tuple,np.ndarray],
                longitude: tp.Union[list,tuple,np.ndarray],
                xlaspl: np.ndarray,
                xlospl: np.ndarray,
                xraspl: np.ndarray):
    """Evaluate spherical splines at a set of locations.

    Parameters
    ----------
    latitude : tp.Union[list,tuple,np.ndarray]
        Latitudes of locations queried
    longitude : tp.Union[list,tuple,np.ndarray]
        Longitudes of locations queried
    xlaspl : np.ndarray
        Latitude locations of splines
    xlospl : np.ndarray
        Longitude locations of splines
    xraspl : np.ndarray
        Radius of splines

    Returns
    -------
    horcof : np.ndarray
        Value of the horizontal coefficents at each location.
        Size of numpy array is [len(latitude) X ncoefhor]

    """

    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)

    if not len(latitude) == len(longitude): raise AssertionError('latitude and longitude should be of same length')
    if not len(xlaspl) == len(xlospl) == len(xraspl): raise AssertionError('xlaspl,xlospl and xraspl should be of same length')

    ncoefhor = len(xlaspl)
    values = sparse.csr_matrix((len(latitude),ncoefhor)) # empty matrix

    for iloc in range(len(latitude)):
        lat = latitude[iloc]
        lon = longitude[iloc]
        #--- make lon go from 0-360
        if lon<0.: lon=lon+360.
        xlospl[np.where(xlospl<0.)]=xlospl[np.where(xlospl<0.)]+360.
        ncon,colind,con = splcon(lat,lon,ncoefhor,xlaspl,xlospl,xraspl)
        rowind = iloc*np.ones(ncon)
        # update values
        values = values + sparse.csr_matrix((con, (rowind, colind)), shape=(len(latitude),ncoefhor))
    return values

@jit(nopython=True)
def splcon(lat:float,lon:float,ncoefhor:int,xlaspl:np.ndarray,xlospl:np.ndarray,xraspl:np.ndarray):
    """Evaluate spherical splines at a given location. Modified from the Fortran code
    splcon.f that is based on Wang and Dahlen (1995) :cite:`Wang:1995fn`.

    Parameters
    ----------
    lat : float
        latitude query
    lon : float
        longitude query
    ncoefhor : int
        number of splines
    xlaspl : np.ndarray
        spline latitudes
    xlospl : np.ndarray
        spline longitudes
    xraspl : np.ndarray
        spline radii

    Returns
    -------
    ncon,icon,con
        number of splines, spline index, and value at each location

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    ncon=0
    con=np.zeros(ncoefhor)
    icon=np.zeros(ncoefhor,dtype=int64)
    for iver in range(ncoefhor):
        if lat>xlaspl[iver]-2.*xraspl[iver]:
            if lat<xlaspl[iver]+2.*xraspl[iver]:
                dd=sind(xlaspl[iver])*sind(lat)
                dd=dd+cosd(xlaspl[iver])*cosd(lat)*cosd(lon-xlospl[iver])
                dd=acosd(dd)
                if dd <= xraspl[iver]*2.:
                    icon[ncon] = iver
                    rn=dd/xraspl[iver]
                    dr=rn-1.
                    if rn <= 1.:
                        con[ncon] = (0.75*rn-1.5)*(rn**2)+1.
                    elif rn > 1.:
                        con[ncon] = ((-0.25*dr+0.75)*dr-0.75)*dr+0.25
                    else:
                        con[ncon] = 0.
                    ncon=ncon+1
    return ncon,icon[:ncon],con[:ncon]

def eval_ylm(latitude: tp.Union[list,tuple,np.ndarray],
             longitude: tp.Union[list,tuple,np.ndarray],
             lmaxhor: int,
             weights: tp.Union[None,list,tuple,np.ndarray] = None,
             grid: bool = False,
             norm: str = 'ylm'):
    """Evaluate spherical harmonics with a specific normalization.

    Parameters
    ----------
    latitude : tp.Union[list,tuple,np.ndarray]
        Latitudes of locations queried
    longitude : tp.Union[list,tuple,np.ndarray]
        Longitudes of locations queried
    lmaxhor : int
        Maximum spherical harmonic degree
    weights : tp.Union[None,list,tuple,np.ndarray], optional
        Weights to multiply the bases with i.e. coefficients, by default None
    grid : bool, optional
        Create a grid based on a combination of latitudes and longitudes, by default False
    norm : str, optional
        Spherical harmonic normalization, by default 'ylm'

    Returns
    -------
    horcof
        Value of the horizontal coefficents at each location.
        Size of numpy array is [len(latitude) X ((lmaxhor+1)^2)]

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # convert to numpy arrays
    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)
    nlat = len(latitude)
    nlon = len(longitude)

    #checks
    if grid:
        nrows,latitude,longitude = makegrid(latitude,longitude)
    else:
        if not (latitude.ndim == longitude.ndim == 1): raise ValueError("latitude, longitude or depth_in_km should be one-dimensional arrays")
        nrows = nlat

    ncoefhor = np.power(lmaxhor+1,2) # numpye of coefficients upto Lmax
    if np.any(weights == None):
        horcof = sparse.csr_matrix((nrows,ncoefhor)) # empty matrix
    else:
        weights = convert2nparray(weights)
        assert(len(weights)==ncoefhor),' length of weights and lat/lon need to be same'
        values = np.zeros_like(latitude)

    if not (len(latitude) == len(longitude)): raise ValueError("latitude, longitude or depth_in_km should be of same length if not making grid = False")
    for iloc in range(len(latitude)):
        lat = latitude[iloc]
        lon = longitude[iloc]
        #--- make lon go from 0-360
        if lon<0.: lon=lon+360.
        # wk1,wk2,wk3 are legendre polynomials of size Lmax+1
        # ylmcof is the value of Ylm
        if norm == 'ylm':
            ylmcof, _ , _ , _ = ylm(lat,lon,lmaxhor,ncoefhor,lmaxhor+1)
        elif norm == 'shold':
            ylmcof = shold(lat,lon,lmaxhor,ncoefhor)
        if np.any(weights == None):
            rowind = iloc*np.ones(ncoefhor)
            colind = np.arange(ncoefhor)
            # update values
            horcof = horcof + sparse.csr_matrix((ylmcof, (rowind, colind)), shape=(nrows,ncoefhor))
        else:
            values[iloc] = np.sum(ylmcof*weights)

    return horcof if np.any(weights == None) else values

def eval_pixel(latitude: tp.Union[list,tuple,np.ndarray],
               longitude: tp.Union[list,tuple,np.ndarray],
               xlapix: tp.Union[list,tuple,np.ndarray],
               xlopix: tp.Union[list,tuple,np.ndarray],
               xsipix: tp.Union[list,tuple,np.ndarray]):
    """Evaluate pixel values at a set of locations.

    Parameters
    ----------
    latitude : tp.Union[list,tuple,np.ndarray]
        Latitudes of locations queried
    longitude : tp.Union[list,tuple,np.ndarray]
        Longitudes of locations queried
    xlapix : tp.Union[list,tuple,np.ndarray]
        Pixel center latitudes
    xlopix : tp.Union[list,tuple,np.ndarray]
        Pixel center longitudes
    xsipix : tp.Union[list,tuple,np.ndarray]
        Pixel sizes in degrees

    Returns
    -------
    horcof
        Value of the horizontal coefficents at each location.
        Size of numpy array is [len(latitude) X len(xsipix)]

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)

    if not len(latitude) == len(longitude): raise AssertionError('latitude and longitude should be of same length')
    if not len(xlapix) == len(xlopix) == len(xsipix): raise AssertionError('xlapix,xlopix,xsipix should be of same length')
    if len(np.unique(xsipix)) > 1:
        strout=''
        for ii in range(len(np.unique(xsipix))): strout = strout+', '+str(np.unique(xsipix)[ii])
        print('Warning: multiple pixel sizes in the PIXEL basis fir evaluation in bases.eval_pixel. Sizes are ')

    labo = xlapix-xsipix/2.; lato = xlapix+xsipix/2.
    labo = np.clip(labo,-90.,90.); lato = np.clip(lato,-90.,90.)
    # go from (-180,180) to (0,360)
    xlopix[np.where(xlopix<0.)]=xlopix[np.where(xlopix<0.)]+360.
    lole = xlopix-xsipix/2.; lori = xlopix+xsipix/2.
    lori = np.clip(lori,0.,360.); lole = np.clip(lole,0.,360.)

    horcof = sparse.csr_matrix((len(latitude),len(xsipix))) # empty matrix
    for iloc,lat in enumerate(latitude):
        lon = longitude[iloc]
        #--- make lon go from 0-360
        if lon<0.: lon=lon+360.
        # check if the location lies within pixel
        if lat != 90.:
            latfind = np.intersect1d(np.where(lat<lato), np.where(lat>=labo))
        else:
            latfind = np.intersect1d(np.where(lat<=lato), np.where(lat>=labo))
        if lon != 360.:
            lonfind = np.intersect1d(np.where(lon<lori), np.where(lon>=lole))
        else:
            lonfind = np.intersect1d(np.where(lon<=lori), np.where(lon>=lole))
        findindex = np.intersect1d(latfind,lonfind)
        if len(findindex) != 1: raise ValueError('found '+str(len(findindex))+' pixels for the location ('+str(lat)+','+str(lon)+')')

        rowind = iloc*np.ones_like(findindex)
        values = np.ones_like(findindex,dtype=float)
        colind = np.array(findindex)
        # update values
        horcof = horcof + sparse.csr_matrix((values, (rowind, colind)), shape=(len(latitude),len(xsipix)))

    return horcof
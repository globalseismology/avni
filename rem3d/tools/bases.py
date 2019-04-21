#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import numpy as np
from collections import Counter
from scipy import sparse

####################### IMPORT REM3D LIBRARIES  #######################################
from rem3d.f2py import vbspl,dbsplrem,ylm
from .trigd import sind,cosd,acosd
from .common import convert2nparray
#######################################################################################

def eval_vbspl(depths,knots):
    """
    Evaluate the cubic spline know with second derivative as 0 at end points.

    Input parameters:
    ----------------

    depth: value or array of depths queried in km

    knots: numpy array or list of depths of spline knots

    Output:
    ------

    vercof, dvercof: value of the spline coefficients at each depth and its derivative.
                      Both arrays have size (Ndepth, Nknots).
    """
    if isinstance(knots, (list,tuple,np.ndarray)):
        knots = np.asarray(knots)
        knots = np.sort(knots)
    else:
        raise TypeError('knots must be list or tuple, not %s' % type(knots))

    # convert to numpy arrays
    depths = convert2nparray(depths)

    # find repeated values
    repeats = [item for item, count in Counter(knots).items() if count > 1]
    repeats_gt_2= [item for item, count in Counter(knots).items() if count > 2]
    if len(repeats_gt_2) != 0: raise ValueError('Cannot have more than 2 repetitions in knots')

    if len(repeats) > 0: # if there are repeated knots, splits it
        split_at = []
        for index,value in enumerate(repeats):
            split_at.append(np.where(knots==value)[0][1])
        knots_list = np.split(knots, split_at)
        for knots in knots_list:
            if len(knots) < 4:
                raise ValueError('Atleast 4 knots need to be defined at or between '+str(min(knots))+' and '+str(max(knots))+' km')

        jj = 0
        for depth in depths:
            jj = jj + 1
            for index,value in enumerate(knots_list):
                # create the arrays as Fortran-contiguous
                splpts = np.array(value.tolist(), order = 'F')
                #Undefined if depth does not lie within the depth extents of knot points
                if depth < min(value) or depth > max(value):
                    temp1 = temp2 = np.zeros_like(splpts)
                else:
                    (temp1, temp2) = vbspl(depth,len(splpts),splpts)
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
    else:
        if len(knots) < 4:
            raise ValueError('Atleast 4 knots need to be defined at or between '+str(min(knots))+' and '+str(max(knots))+' km')
        # create the arrays as Fortran-contiguous
        splpts = np.array(knots.tolist(), order = 'F')
        jj = 0
        for depth in depths:
            jj = jj + 1
            #Undefined if depth does not lie within the depth extents of knot points
            if depth < min(knots) or depth > max(knots):
                vercof_temp = dvercof_temp = np.zeros_like(splpts)
            else:
                (vercof_temp, dvercof_temp) = vbspl(depth,len(splpts),splpts)
            if jj == 1:
                vercof = vercof_temp; dvercof = dvercof_temp
            else:
                vercof = np.vstack([vercof,vercof_temp])
                dvercof = np.vstack([dvercof,dvercof_temp])
    return vercof, dvercof


def eval_splrem(radius, radius_range, nsplines):
    """
    Evaluate the cubic spline know with second derivative as 0 at end points.

    Input parameters:
    ----------------

    radius: value or array of radii queried

    radius_range: limits of the radius limits of the region

    nsplines: number of splines within the range

    Output:
    ------

    vercof, dvercof: value of the polynomial coefficients at each depth and derivative.
                      Both arrays have size (Nradius, Nsplines).
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


def eval_polynomial(radius, radius_range, rnorm, types = None):
    """
    Evaluate the cubic spline know with second derivative as 0 at end points.

    Input parameters:
    ----------------

    radius: value or array of radii queried

    radius_range: limits of the radius limits of the region

    types: polynomial coefficients to be used for calculation. Options are : TOP,
                  TOP, BOTTOM, CONSTANT, LINEAR, QUADRATIC, CUBIC. 
                  default: ['CONSTANT','LINEAR']

    rnorm: normalization for radius, usually the radius of the planet

    Output:
    ------

    vercof : value of the polynomial coefficients at each depth, size (Nradius).
    """
    #defaults
    if types is None: types= ['CONSTANT','LINEAR']

    # convert to numpy arrays
    radiusin = convert2nparray(radius)
    radius_range = convert2nparray(radius_range)

    if radius_range.shape[1] != 2 or not isinstance(radius_range, (list,tuple,np.ndarray)):
        raise TypeError('radius_range must be list , not %s' % type(radius_range))

    # keys in coefficients should be acceptable
    choices = ['TOP', 'BOTTOM', 'CONSTANT', 'LINEAR', 'QUADRATIC', 'CUBIC']
    if not np.all([key in choices for key in types]): raise AssertionError()
    npoly = len(radius_range)*len(types)
    # first find whether CONSTANT/LINEAR or TOP/BOTTOM
    for ii in range(len(radius_range)):
        if not np.all(np.sort(radius_range[ii])==radius_range[ii]): raise AssertionError('radius_range needs to be sorted')
    rbot=radius_range[0]/rnorm
    rtop=radius_range[1]/rnorm
    findtopbot = np.any([key in ['BOTTOM','TOP'] for key in types])
    findconstantlinear = np.any([key in ['CONSTANT','LINEAR'] for key in types])

    if findtopbot and findconstantlinear: raise ValueError('Cannot have both BOTTOM/TOP and CONSTANT/LINEAR as types in eval_polynomial')

    for irad in range(len(radiusin)):
        temp = np.zeros(npoly)
        dtemp = np.zeros(npoly)
        for irange in range(len(radius_range)):
            #Undefined if depth does not lie within the depth extents of knot points
            if radiusin[irad] <= min(radius_range[irange]) or radiusin[irad] > max(radius_range[irange]):
                # <= so that the boundary depth belongs to only one radial kernel
                for itype in range(len(types)):
                    ii = irange*len(types)+itype
                    temp[ii]=0.
                    dtemp[ii]=0.
            else:
                rn=radiusin[irad]/rnorm
                for itype in range(len(types)):
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

def eval_splcon(latitude,longitude,xlaspl,xlospl,xraspl):
    """
    Evaluate the continuous lateral splines.

    Input parameters:
    ----------------

    latitude,longitude: location queried

    xlaspl,xlospl,xraspl: loocation and radius of splines

    Output:
    ------

    horcof : value of the horizontal coefficents at each location.
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
        ncon,icon,con = splcon(lat,lon,ncoefhor,xlaspl,xlospl,xraspl)
        rowind = iloc*np.ones(ncon)
        colind = []
        for ii in range(ncon): colind.append(icon[ii])
        colind = np.array(colind)
        # update values
        values = values + sparse.csr_matrix((con[:ncon], (rowind, colind)), shape=(len(latitude),ncoefhor))
    return values

def splcon(lat,lon,ncoefhor,xlaspl,xlospl,xraspl):
    ncon=0;con=[];icon=[]
    for iver in range(ncoefhor):
        if lat>xlaspl[iver]-2.*xraspl[iver]:
            if lat<xlaspl[iver]+2.*xraspl[iver]:
                dd=sind(xlaspl[iver])*sind(lat)
                dd=dd+cosd(xlaspl[iver])*cosd(lat)*cosd(lon-xlospl[iver])
                dd=acosd(dd)
                if dd <= xraspl[iver]*2.:
                    ncon=ncon+1
                    icon.append(iver)
                    rn=dd/xraspl[iver]
                    dr=rn-1.
                    if rn <= 1.:
                        con.append((0.75*rn-1.5)*(rn**2)+1.)
                    elif rn > 1.:
                        con.append(((-0.25*dr+0.75)*dr-0.75)*dr+0.25)
                    else:
                        con.append(0.)
    con=np.array(con);icon=np.array(icon)
    return ncon,icon,con

def eval_ylm(latitude,longitude,lmaxhor):
    """
    Evaluate spherical harmonics.

    Input parameters:
    ----------------

    latitude,longitude: location queried

    lmaxhor: maximum spherical harmonic degree

    Output:
    ------

    horcof : value of the horizontal coefficents at each location.
             Size of numpy array is [len(latitude) X ((lmaxhor+1)^2)]
    """

    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)

    if not len(latitude) == len(longitude): raise AssertionError('latitude and longitude should be of same length')
    ncoefhor = np.power(lmaxhor+1,2) # numpye of coefficients upto Lmax
    horcof = sparse.csr_matrix((len(latitude),ncoefhor)) # empty matrix
    for iloc in range(len(latitude)):
        lat = latitude[iloc]
        lon = longitude[iloc]
        #--- make lon go from 0-360
        if lon<0.: lon=lon+360.
        # wk1,wk2,wk3 are legendre polynomials of size Lmax+1
        # ylmcof is the value of Ylm
        ylmcof, _ , _ , _ = ylm(lat,lon,lmaxhor,ncoefhor,lmaxhor+1)
        rowind = iloc*np.ones(ncoefhor)
        colind = np.arange(ncoefhor)
        # update values
        horcof = horcof + sparse.csr_matrix((ylmcof, (rowind, colind)), shape=(len(latitude),ncoefhor))
    return horcof

def eval_pixel(latitude,longitude,xlapix,xlopix,xsipix):
    """
    Evaluate spherical harmonics.

    Input parameters:
    ----------------

    latitude,longitude: location queried

    xlapix,xlopix,xsipix: loocation and size of pixels

    Output:
    ------

    horcof : value of the horizontal coefficents at each location.
             Size of numpy array is [len(latitude) X len(xsipix)]
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
    lole = xlopix-xsipix/2.; lori = xlopix+xsipix/2.
    lori[np.where(lori<0.)]=lori[np.where(lori<0.)]+360.
    lole[np.where(lole<0.)]=lole[np.where(lole<0.)]+360.
    lori[np.where(lori>360.)]=lori[np.where(lori>360.)]-360.
    lole[np.where(lole>360.)]=lole[np.where(lole>360.)]-360.
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

        rowind = iloc*np.ones_like(findindex)
        values = np.ones_like(findindex,dtype=np.float)
        colind = np.array(findindex)
        # update values
        horcof = horcof + sparse.csr_matrix((values, (rowind, colind)), shape=(len(latitude),len(xsipix)))

    return horcof



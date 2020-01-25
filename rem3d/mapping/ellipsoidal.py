#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function

import multiprocessing
from joblib import Parallel, delayed
import numpy as np
import pdb
############################### PLOTTING ROUTINES ################################
from ..tools.trigd import atand,tand
from ..tools.common import convert2nparray
from ..f2py import ddelazgc # geolib library from NSW
from .. import constants
###############################

def get_distaz(eplat,eplon,stlat,stlon,num_cores=1):
    """Get the distance and azimuths from positions in geographic coordinates"""

    geoco = constants.geoco.magnitude
    eplat = convert2nparray(eplat)
    eplon = convert2nparray(eplon)
    stlat = convert2nparray(stlat)
    stlon = convert2nparray(stlon)
    if not (len(stlat) == len(stlon) == len(eplat) == len(eplon)):
        pdb.set_trace()
        raise ValueError('latitude and longitude need to be of same length')

    if len(eplat) > 1: # if the input is a list loop
        delta=[];azep=[];azst=[]
        # Standard checks on number of cores
        avail_cores = multiprocessing.cpu_count()
        if num_cores > avail_cores:
            raise ValueError("Number of cores requested ("+str(num_cores)+") is higher than available ("+str(avail_cores)+")")
        # Crate list of tuples of job arguments and pass to parallel using a helper routine
        job_args = [(atand(geoco*tand(item_lat)),eplon[jj],atand(geoco*tand(stlat[jj])),stlon[jj]) for jj, item_lat in enumerate(eplat)]
        temp=Parallel(n_jobs=num_cores)(delayed(delazgc_helper)(ii) for ii in job_args)
        for il in temp: delta.append(il[0]);azep.append(il[1]);azst.append(il[2])

    else:
        elat=atand(geoco*tand(eplat))
        elon=eplon
        slat=atand(geoco*tand(stlat))
        slon=stlon
        delta,azep,azst = ddelazgc(elat,elon,slat,slon)

    return delta,azep,azst

def delazgc_helper(args):
    return ddelazgc(*args)

def geographic_to_geocentric(latin):
    """
    convert a geographic coordinate to geocentric coordinate
    """
    xlat=latin
    fac = constants.geoco.magnitude
    theta = radians(90.0-xlat)
    theta = pi/2.-atan2(fac*cos(theta),sin(theta))
    latout=90.0-degrees(theta)

    return latout

def geocentric_to_geographic(latin):
    """
    convert a geocentric coordinate to geographic coordinate
    """
    xlat=latin
    fac = constants.geoco.magnitude
    if xlat != 0.:
        theta = radians(90.0-xlat)
        theta = atan(fac/tan(pi/2.-theta))
        latout=90.0-degrees(theta)
    else:
        latout=xlat
    if xlat < 0.: latout = latout - 180.

    return latout

def inpolygon(latitude,longitude,polygon_latitude,polygon_longitude,num_cores=1, orientation = 'anti-clockwise', threshold = 0.01):
    """
    Finds whether a (set of) point(s) is within a closed polygon

    Input Parameters:
    ----------------

    latitude,longitude: set of queried locations

    polygon_latitude,polygon_longitude: closed points that define the polygon.
                                First and last points need to be the same.

    num_cores: Number of cores to use for the calculations

    orientation: clockwise or anti-clockwise orientation of points specified above

    threshold: limit to which the sum of azimuth check to (-)360 degrees is permitted
               to be defined as within the polygon.

    Output:
    ------

    within: logical array containing the same number of elements as latitude/longitude
    """
    # convert to numpy arrays. Various checks
    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)
    if len(latitude) != len(longitude): raise ValueError('latitude and longitude need to be of same length')
    polylat = convert2nparray(polygon_latitude)
    polylon = convert2nparray(polygon_longitude)
    if len(polylat) != len(polylon): raise ValueError('polygon latitude and longitude need to be of same length')
    if (polygon_latitude[-1] != polygon_latitude[0]) or (polygon_longitude[-1] != polygon_longitude[0]): raise ValueError('polygon should be closed and therefore have same start and end points')
    nvertices = len(polylon)
    if orientation not in ['clockwise','anti-clockwise']: raise ValueError('orientation needs to be clockwise or anti-clockwise')

    # define the output array; assume everything is outside
    within = np.zeros_like(longitude,dtype=bool)

    # loop over all queries points
    for ii,lat in enumerate(latitude):
        lon = longitude[ii]
        # Get the azimuth of this location to all vertices of the polygon
        dist,azep,azst = get_distaz(np.repeat(lat,nvertices),np.repeat(lon,nvertices),polylat,polylon,num_cores=num_cores)

    #   Since the polygons are anti-clockwise, check for -360 instead of 360
    #   360 gives the plate for the anitpodal location in this case
        angle=0.
        for jj, azimuth in enumerate(azep):
            if jj > 0:
                add = azimuth - azold
                if add > 180.: add=add-360.
                if add < -180.: add=add+360.
                angle=angle+add
            azold = azimuth
        if orientation == 'anti-clockwise':
            if abs(angle+360) < threshold: within[ii] = True
        elif orientation == 'clockwise':
            if abs(angle-360) < threshold: within[ii] = True
    return within[0] if len(within)==1 else within

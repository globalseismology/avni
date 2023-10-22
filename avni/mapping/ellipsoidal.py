#!/usr/bin/env python

#####################  IMPORT STANDARD MODULES   #########################

# python 3 compatibility
from __future__ import absolute_import, division, print_function

import multiprocessing
from joblib import Parallel, delayed
import numpy as np
import typing as tp
from math import cos, sin, tan, atan, atan2, degrees, radians

####################### IMPORT AVNI LIBRARIES  ###########################

from ..tools.trigd import atand,tand
from ..tools.common import convert2nparray
from ..f2py import ddelazgc # geolib library from NSW
from .. import constants

##########################################################################

def get_distaz(eplat: tp.Union[float,list,tuple,np.ndarray],
               eplon: tp.Union[float,list,tuple,np.ndarray],
               stlat: tp.Union[float,list,tuple,np.ndarray],
               stlon: tp.Union[float,list,tuple,np.ndarray],
               num_cores: int = 1) -> (tp.Union[np.ndarray,float],
               tp.Union[np.ndarray,float],tp.Union[np.ndarray,float]):
    """Get the distance and azimuths between pairs of positions in geographic coordinates

    This function first converts the queried station and source locations to
    geocentric coordinates. In practice, we projects all points from an ellipsoid
    of flatness (f) to a sphere of equatorial radius (a_e) though the geocentric
    conversion factor (W or geoco below). This is also called the parametric or reduced
    latitude conversion, introduced by Legendre and Bessel who solved problems
    for geodesics on the ellipsoid by transforming them to an equivalent problem for
    spherical geodesics by using this smaller geocentric latitude.

    The spherical law of cosines formula is then used to calculate distances on this
    spherical geodesic of radius a_e. This procedure stretches the angular distances
    between adjacent geographic latitudes nearer to the poles and equator, which imitates the
    behavior in an ellipsoid where adjacent latitudes become finely spaced.
    The equatorial radius is used instead of mean radius (R) for conversion of distances
    to km since this conversion is strictly valid for a sphere of radius a_e
    and in order to obtain accurate distances near the equator.

    Parameters
    ----------
    eplat : tp.Union[float,list,tuple,np.ndarray]
        Latitudes of source location(s) in Geographic Coordinates
    eplon : tp.Union[float,list,tuple,np.ndarray]
        Longitudes of source location(s) in Geographic Coordinates
    stlat : tp.Union[float,list,tuple,np.ndarray]
        Latitudes of station location(s) in Geographic Coordinates
    stlon : tp.Union[float,list,tuple,np.ndarray]
        Longitudes of station location(s) in Geographic Coordinates
    num_cores : int, optional
        Number of cores to use in the calculation, by default 1

    Returns
    -------
    delta,azep,azst
        A tuple with elements as distance in degrees, azimuth from source(s), and
        backazimuth from station(s).

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    geoco = constants.geoco.magnitude
    eplat = convert2nparray(eplat)
    eplon = convert2nparray(eplon)
    stlat = convert2nparray(stlat)
    stlon = convert2nparray(stlon)
    if not (len(stlat) == len(stlon) == len(eplat) == len(eplon)):
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
    """A helper function to parallelize ddelazgc"""
    return ddelazgc(*args)

def geographic_to_geocentric(latin: float) -> float:
    """Convert a geographic latitude to geocentric latitude

    Parameters
    ----------
    latin : float
        Input latitude in geographic coordinate

    Returns
    -------
    float
        Output geocentric latitude

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    xlat=latin
    fac = constants.geoco.magnitude
    theta = radians(90.0-xlat)
    theta = np.pi/2.-atan2(fac*cos(theta),sin(theta))
    latout=90.0-degrees(theta)

    return latout

def geocentric_to_geographic(latin: float) -> float:
    """Convert a geocentric coordinate to geographic coordinate

    Parameters
    ----------
    latin : float
        Input latitude in geocentric coordinate

    Returns
    -------
    float
        Output geographic latitude

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    xlat=latin
    fac = constants.geoco.magnitude
    if xlat != 0.:
        theta = radians(90.0-xlat)
        theta = atan(fac/tan(np.pi/2.-theta))
        latout=90.0-degrees(theta)
    else:
        latout=xlat
    if xlat < 0.: latout = latout - 180.

    return latout

def inpolygon(latitude: tp.Union[float,list,tuple,np.ndarray],
              longitude: tp.Union[float,list,tuple,np.ndarray],
              polygon_latitude: tp.Union[list,tuple,np.ndarray],
              polygon_longitude: tp.Union[list,tuple,np.ndarray],
              num_cores: int = 1, orientation: str = 'anti-clockwise',
              threshold: float = 1E-6) -> tp.Union[np.ndarray,bool]:
    """Finds whether a (set of) point(s) is(are) within a closed polygon

    Parameters
    ----------
    latitude,longitude : tp.Union[float,list,tuple,np.ndarray]
        Set of queried locations in Geographic Coordinates
    polygon_latitude,polygon_longitude : tp.Union[list,tuple,np.ndarray]
        Closed points that define the polygon.
        First and last points need to be the same.
    num_cores : int, optional
        Number of cores to use for the calculations, by default 1
    orientation : str, optional
        clockwise or anti-clockwise orientation of points specified above, by default 'anti-clockwise'
    threshold : float, optional
        Limit to which the sum of azimuth check to (-)360 degrees is permitted
        to be defined as within the polygon, by default 1E-6

    Returns
    -------
    tp.Union[np.ndarray,bool]
        Logical array of bool(s) containing the same number of elements as
        latitude/longitude

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
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
    #   Also 0 means the same as 360/-360 so that is also checked
        angle=0.
        for jj, azimuth in enumerate(azep):
            if jj > 0:
                add = azimuth - azold
                if add > 180.: add=add-360.
                if add < -180.: add=add+360.
                angle=angle+add
            azold = azimuth
        if orientation == 'anti-clockwise':
            if abs(angle+360) <= threshold or abs(angle) <= threshold: within[ii] = True
        elif orientation == 'clockwise':
            if abs(angle-360) <= threshold or abs(angle) <= threshold: within[ii] = True
    return within[0] if len(within)==1 else within
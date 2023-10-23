#!/usr/bin/env python

#####################  IMPORT STANDARD MODULES   #########################

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int

import numpy as np #for numerical analysis
import typing as tp

####################### IMPORT AVNI LIBRARIES  ###########################

from ..tools.common import convert2nparray
from .. import constants

##########################################################################

def intersection(path1start: tp.Union[list,np.ndarray],
                 path1brngEnd: tp.Union[float,int,list,np.ndarray],
                 path2start: tp.Union[list,np.ndarray],
                 path2brngEnd: tp.Union[float,int,list,np.ndarray]):
    """Get the intersection of two great circle paths. Input can either be
    point and bearings or two sets of points.

    If c1 & c2 are great circles through start and end points
    (or defined by start point + bearing),
    then candidate intersections are simply c1 X c2 & c2 x c1
    most of the work is deciding correct intersection point to select!
    If bearing is given, that determines which intersection,
    if both paths are defined by start/end points, take closer intersection
    https://www.movable-type.co.uk/scripts/latlong-vectors.html#intersection

    Parameters
    ----------
    path1start : tp.Union[list,np.ndarray]
        Location of start point of the first curve in [latitude, longitude]
    path1brngEnd : tp.Union[float,int,list,np.ndarray]
        End point of first curve either in terms of a bearing (float/int) or location [latitude, longitude]
    path2start : tp.Union[float,int,list,np.ndarray]
        Location of start point of the second curve in [latitude, longitude]
    path2brngEnd : tp.Union[float,int,list,np.ndarray]
        End point of second curve either in terms of a bearing (float/int) or location [latitude, longitude]

    Returns
    -------
    intersection
        Preferred intersection of the two curves based on the input

    antipode
        Antipode of the intersection where the great-circle curves meet as well

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # find out what the types are
    # c1 & c2 are vectors defining great circles through start & end points
    # p X c gives initial bearing vector
    if isinstance(path1brngEnd, (float,int)):
        path1def = 'bearing' # path 1 defined by endpoint
        p1end = spher2cart(getDestination(path1start[0],path1start[1],path1brngEnd,180.))
    else:
        path1def = 'endpoint' #path 1 defined by initial bearing
        p1end = path1brngEnd
    if isinstance(path2brngEnd, (float,int)):
        path2def = 'bearing'
        p2end = spher2cart(getDestination(path2start[0],path2start[1],path2brngEnd,180.))
    else:
        path2def = 'endpoint'
        p2end = path2brngEnd
    case = path1def + '+' + path2def

    # convert to spherical coordinates
    p1 = spher2cart(path1start)
    p2 = spher2cart(path2start)

    # Get normal to planes containing great circles
    # np.cross product of vector to each point from the origin
    N1 = np.cross(p1, p1end)
    N2 = np.cross(p2, p2end)

    # Find line of intersection between two planes
    L = np.cross(N1, N2)

    # Find two intersection points
    X1 = L / np.sqrt(L[0]**2 + L[1]**2 + L[2]**2)
    X2 = -X1

    #Check if correct
    if np.any(np.isnan(X1)): raise ValueError('not found an intersection point between ',path1start,' with azimuth ',path1brngEnd,' and ', path2start, ' with azimuth ',path2brngEnd)

    # convert back to spherical
    i1 = cart2spher(X1)[1:]
    i2 = cart2spher(X2)[1:]

    # selection of intersection point depends on
    # how paths are defined (bearings or endpoints)
    if case == 'bearing+bearing':
        #if NXp1.i1 is +ve, the initial bearing is towards i1
        # otherwise towards antipodal i2
        dir1 = np.sign(np.dot(np.cross(N1,p1),X1)) #c1Xp1.X1 +ve means p1 bearing points to X1
        dir2 = np.sign(np.dot(np.cross(N2,p2),X1)) #c2Xp2.X1 +ve means p2 bearing points to X1
        if dir1 + dir2 == 2: # dir1, dir2 both +ve, 1 & 2 both pointing to X1
            intersect = i1
            antipode = i2
        elif dir1 + dir2 == -2: #dir1, dir2 both -ve, 1 & 2 both pointing to X2
            intersect = i2
            antipode = i1
        elif dir1 + dir2 == 0:
            # dir1, dir2 opposite; intersection is at further-away intersection point
            # take opposite intersection from mid-point of p1 & p2 [is this always true?]
            if np.dot(p1+p2,X1) > 0 :
                intersect = i2
                antipode = i1
            else:
                intersect = i1
                antipode = i2
    elif  case == 'bearing+endpoint':  #use bearing c1 X p1
        dir1 = np.sign(np.dot(np.cross(N1,p1),X1)) #c1Xp1.X1 +ve means p1 bearing points to X1
        if dir1 > 0:
            intersect = i1
            antipode = i2
        else:
            intersect = i2
            antipode = i1
    elif  case == 'endpoint+bearing': #use bearing c2 X p2
        dir2 = np.sign(np.dot(np.cross(N2,p2),X1)) #c2Xp2.X1 +ve means p2 bearing points to X1
        if dir2 > 0:
            intersect = i1
            antipode = i2
        else:
            intersect = i2
            antipode = i1
    elif case == 'endpoint+endpoint': #select nearest intersection to mid-point of all points
        mid = p1+p2+p1end+p2end
        if np.dot(mid,X1) > 0 :
            intersect = i1
            antipode = i2
        else:
            intersect = i2
            antipode = i1

    ## Attempted with pygeodesy. Has issues for certain configurations.
#     if path1start[1] > 180.: path1start[1] = path1start[1] -360.
#     path1start = LatLon(path1start[0],path1start[1])
#     if path2start[1] > 180.: path2start[1] = path2start[1] -360.
#     path2start = LatLon(path2start[0],path2start[1])
#     # find out what the types are
#     # c1 & c2 are vectors defining great circles through start & end points
#     # p X c gives initial bearing vector
#     if isinstance(path1brngEnd, (float,int)):
#         path1def = 'bearing' # path 1 defined by endpoint
#     else:
#         path1def = 'endpoint' #path 1 defined by initial bearing
#         if path1brngEnd[1] > 180.: path1brngEnd[1] = path1brngEnd[1] -360.
#         path1brngEnd = LatLon(path1brngEnd[0],path1brngEnd[1])
#     if isinstance(path2brngEnd, (float,int)):
#         path2def = 'bearing'
#     else:
#         path2def = 'endpoint'
#         if path2brngEnd[1] > 180.: path2brngEnd[1] = path2brngEnd[1] -360.
#         path2brngEnd = LatLon(path2brngEnd[0],path2brngEnd[1])
#     intersect = path1start.intersection(path1brngEnd, path2start, path2brngEnd)

    return intersect,antipode

def midpoint(lat1: tp.Union[int,float], lon1: tp.Union[int,float],
             lat2: tp.Union[int,float], lon2: tp.Union[int,float]):
    """Get the mid-point from positions in geographic coordinates. Input values as degrees"""
    from pygeodesy.sphericalNvector import LatLon

    if lon1 > 180.: lon1 = lon1 - 360.
    if lon2 > 180.: lon2 = lon2 - 360.

    start = LatLon(lat1,lon1)
    end = LatLon(lat2,lon2)
    mid = start.midpointTo(end)
    return [mid.lat, mid.lon]

def cart2spher(xyz: tp.Union[list,np.ndarray]) -> np.ndarray:
    """Convert from cartesian to spherical coordinates.

    Parameters
    ----------
    xyz : tp.Union[list,np.ndarray]
        Data containing columns for x,y,z locations in cartesian coordinates

    Returns
    -------
    np.ndarray
        Output containing radius, latitude and longitude in spherical coordinates

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    xyz = convert2nparray(xyz)
    rlatlon = np.zeros(xyz.shape)
    if xyz.ndim == 2:
        xy = xyz[:,0]**2 + xyz[:,1]**2
        rlatlon[:,0] = np.sqrt(xy + xyz[:,2]**2)
    #     ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
        #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
        rlatlon[:,1] = np.arctan2(xyz[:,2], np.sqrt(xy))/np.pi*180. # latitude
        rlatlon[:,2] = np.arctan2(xyz[:,1], xyz[:,0])/np.pi*180.
    elif xyz.ndim == 1:
        xy = xyz[0]**2 + xyz[1]**2
        rlatlon[0] = np.sqrt(xy + xyz[2]**2)
        rlatlon[1] = np.arctan2(xyz[2], np.sqrt(xy))/np.pi*180. # latitude
        rlatlon[2] = np.arctan2(xyz[1], xyz[0])/np.pi*180.
    else:
        raise ValueError('dimension of xyz should be 1 or 2')
    return rlatlon

def spher2cart(rlatlon: tp.Union[list,np.ndarray]) -> np.ndarray:
    """Convert from spherical to cartesian coordinates.

    Parameters
    ----------
    rlatlon : tp.Union[list,np.ndarray]
        Data containing columns for x,y,z locations in spherical coordinates

    Returns
    -------
    np.ndarray
        Output containing x,y,z in spherical coordinates

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    rlatlon = convert2nparray(rlatlon)
    xyz = np.zeros(rlatlon.shape)
    if rlatlon.ndim == 2:
        # assumes unit sphere if r is not provided
        if rlatlon.shape[1] == 2: rlatlon = np.insert(rlatlon,0,1.,axis=1)
        colatitude=90.-rlatlon[:,1]
        xyz[:,0] = rlatlon[:,0]*np.cos(np.pi/180.*rlatlon[:,2])*np.sin(np.pi/180.*colatitude)
        xyz[:,1] = rlatlon[:,0]*np.sin(np.pi/180.*rlatlon[:,2])*np.sin(np.pi/180.*colatitude)
        xyz[:,2] = rlatlon[:,0]*np.cos(np.pi/180.*colatitude)
    elif rlatlon.ndim == 1:
        # assumes unit sphere if r is not provided
        if len(rlatlon) == 2: rlatlon = np.insert(rlatlon,0,1.)
        colatitude=90.-rlatlon[1]
        xyz[0] = rlatlon[0]*np.cos(np.pi/180.*rlatlon[2])*np.sin(np.pi/180.*colatitude)
        xyz[1] = rlatlon[0]*np.sin(np.pi/180.*rlatlon[2])*np.sin(np.pi/180.*colatitude)
        xyz[2] = rlatlon[0]*np.cos(np.pi/180.*colatitude)
    else:
        raise ValueError('dimension of rlatlon should be 1 or 2')
    return xyz

def polar2cart(rtheta: tp.Union[list,np.ndarray]) -> np.ndarray:
    """Convert from polar to cartesian coordinates

    Parameters
    ----------
    rtheta : tp.Union[list,np.ndarray]
        A single r,theta or an array of these for conversion

    Returns
    -------
    np.ndarray
        Output containing x,y in cartesian coordinates

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    rtheta = convert2nparray(rtheta)
    xy = np.zeros(rtheta.shape)
    if rtheta.ndim == 2:
        xy[:,0] = rtheta[:,0]*np.cos(np.pi/180.*rtheta[:,1])
        xy[:,1] = rtheta[:,0]*np.sin(np.pi/180.*rtheta[:,1])
    elif rtheta.ndim == 1:
        xy[0] = rtheta[0]*np.cos(np.pi/180.*rtheta[1])
        xy[1] = rtheta[0]*np.sin(np.pi/180.*rtheta[1])
    else:
        raise ValueError('dimension of rtheta should be 1 or 2')
    return xy

def cart2polar(xy: tp.Union[list,np.ndarray]) -> np.ndarray:
    """Convert from polar to cartesian coordinates

    Parameters
    ----------
    xy : tp.Union[list,np.ndarray]
        A single x,y or an array of these ffor conversion

    Returns
    -------
    np.ndarray
        Output containing r,theta in polar coordinates

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    xy = convert2nparray(xy)
    rtheta = np.zeros(xy.shape)
    if xy.ndim == 2:
        rtheta[:,0] = np.sqrt(np.power(xy[:,1],2)+np.power(xy[:,1],2))
        rtheta[:,1] = np.arctan2(xy[:,1],xy[:,0]) * 180 / np.pi
    elif xy.ndim == 1:
        rtheta[0] = np.sqrt(np.power(xy[1],2)+np.power(xy[1],2))
        rtheta[1] = np.arctan2(xy[1],xy[0]) * 180 / np.pi
    else:
        raise ValueError('dimension of xy should be 1 or 2')
    return rtheta

def getDestination(lat: tp.Union[int,float],lon: tp.Union[int,float],
                   azimuth: tp.Union[int,float], distance) -> list:
    """Returns the location of destination point
    given the start lat, long, aziuth, and distance (in meters).

    Parameters
    ----------
    lat, lon : tp.Union[int,float]
        Start point latitude and longitude
    azimuth : tp.Union[int,float]
        Azimuth to destination
    distance
        Distance to destination (in meters)

    Returns
    -------
    list
        List containing latitude and longitude of destination point
    """
    ''''''
    from pygeodesy.sphericalNvector import LatLon

    R = constants.R.to_base_units().magnitude #Radius of the Earth in m
    if not isinstance(distance, (float,int)): distance = distance.to_base_units().magnitude #Distance m
    if lon > 180.: lon = lon -360.
    start = LatLon(lat,lon)
    end = start.destination(distance,azimuth,R)
    return[end.lat, end.lon]

def calculateBearing(lat1: tp.Union[int,float], lon1: tp.Union[int,float],
                     lat2: tp.Union[int,float], lon2: tp.Union[int,float]):
    """Calculates the azimuth in degrees from start point to end point.

    Parameters
    ----------
    lat1,lon1 : tp.Union[int,float]
        Start point latitude and longitude
    lat2,lon2 : tp.Union[int,float]
        End point latitude and longitude

    Returns
    -------
    bearing
        Bearing to second point

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    from pygeodesy.sphericalNvector import LatLon

    if lon1 > 180.: lon1 = lon1 -360.
    start = LatLon(lat1,lon1)
    if lon2 > 180.: lon2 = lon2 -360.
    end = LatLon(lat2,lon2)
    bearing = start.initialBearingTo(end)
    return bearing

def getIntermediate(lat1: tp.Union[int,float],lon1: tp.Union[int,float],
                    azimuth: tp.Union[int,float],distance: tp.Union[int,float],
                    interval: tp.Union[int,float]) -> list:
    """Returns intermediate coordinates between two coordinate pairs given a desired interval

    Parameters
    ----------
    lat1,lon1 : tp.Union[int,float]
        Start point latitude and longitude
    azimuth : tp.Union[int,float]
        Azimuth to destination
    distance
        Distance to destination (in meters)
    interval : tp.Union[int,float]
        Interval for intermediate points (in meters)

    Returns
    -------
    list
        Coordinates (lat,lon) of intermdiate points
    """
    if lon1 > 180.: lon1 = lon1 -360.
    steps = int(distance / interval)
    coords = []
    coords.append([lat1,lon1])
    for step in np.arange(steps):
        counter = float(interval) * float(step+1)
        coord = getDestination(lat1,lon1,azimuth,counter)
        coords.append(coord)
    return coords

def calculateDistance(lat1: tp.Union[int,float],lon1: tp.Union[int,float],
                      lat2: tp.Union[int,float],lon2: tp.Union[int,float],
                      final_units='m') -> float:
    """Calculates distance between two coordinates

    Parameters
    ----------
    lat1,lon1 : tp.Union[int,float]
        Start point latitude and longitude
    lat2,lon2 : tp.Union[int,float]
        End point latitude and longitude
    final_units : str, optional
        Output distance in km/m/deg, by default 'm'

    Returns
    -------
    float
        Distance between coordinate pairs in m, km or deg
    """
    from pygeodesy.sphericalNvector import LatLon

    if lon1 > 180.: lon1 = lon1 -360.
    if lon2 > 180.: lon2 = lon2 -360.
    start = LatLon(lat1,lon1)
    end = LatLon(lat2,lon2)
    dist = start.distanceTo(end)
    if final_units=='m':
        return dist
    elif final_units=='km':
        gcdist = dist / 1000.
        return gcdist
    elif final_units=='deg':
        gcdist = dist / constants.deg2m
        return gcdist
    else:
        raise ValueError('Final units need to be m, km or deg')
#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import numpy as np #for numerical analysis
from pygeodesy.sphericalNvector import LatLon
import pdb    #for the debugger 
if (sys.version_info[:2] < (3, 0)): range = xrange

############################### PLOTTING ROUTINES ################################        
from ..tools.common import convert2nparray
from .. import constants
###############################
 
def intersection(path1start, path1brngEnd, path2start, path2brngEnd):
    """
    Get the intersection of two great circle paths. Input can either be
    point and bearings or two sets of points
    
    if c1 & c2 are great circles through start and end points 
    (or defined by start point + bearing),
    then candidate intersections are simply c1 X c2 & c2 x c1
    most of the work is deciding correct intersection point to select! 
    if bearing is given, that determines which intersection, 
    if both paths are defined by start/end points, take closer intersection
    https://www.movable-type.co.uk/scripts/latlong-vectors.html#intersection
    
    Input Parameters:
    ----------------
    
    path1start: Location of start point of a curve in [latitude, longitude]
    
    path1brngEnd: End point either in terms of a bearing (float/int) or location
    
    path2start: Similar to path1start for the second curve
    
    path2brngEnd: Similar to path1brngEnd for the second curve
    
    Return:
    ------
    
    intersection: preferred intersection of the two curves based on the input
    
    antipode: antipode of the intersection where the great-circle curves meet as well
    
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
        
def midpoint(lat1, lon1, lat2, lon2):
    """Get the mid-point from positions in geographic coordinates.Input values as degrees"""
    if lon1 > 180.: lon1 = lon1 - 360.
    if lon2 > 180.: lon2 = lon2 - 360.
    start = LatLon(lat1,lon1)
    end = LatLon(lat2,lon2)
    mid = start.midpointTo(end)
    return [mid.lat, mid.lon]
    
def cart2spher(xyz):
    """Convert from cartesian to spherical coordinates
    http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node42.html
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
 
def spher2cart(rlatlon):
    """Convert from spherical to cartesian coordinates
    http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node42.html
    """
    rlatlon = convert2nparray(rlatlon)
    if rlatlon.ndim == 2:
        # assumes unit sphere if r is not provided
        if rlatlon.shape[1] == 2: rlatlon = np.insert(rlatlon,0,1.,axis=1)
        colatitude=90.-rlatlon[:,1]
        xyz = np.zeros(rlatlon.shape)
        xyz[:,0] = rlatlon[:,0]*np.cos(np.pi/180.*rlatlon[:,2])*np.sin(np.pi/180.*colatitude)
        xyz[:,1] = rlatlon[:,0]*np.sin(np.pi/180.*rlatlon[:,2])*np.sin(np.pi/180.*colatitude)
        xyz[:,2] = rlatlon[:,0]*np.cos(np.pi/180.*colatitude)
    elif rlatlon.ndim == 1:
        # assumes unit sphere if r is not provided
        if len(rlatlon) == 2: rlatlon = np.insert(rlatlon,0,1.)
        colatitude=90.-rlatlon[1]
        xyz = np.zeros(rlatlon.shape)
        xyz[0] = rlatlon[0]*np.cos(np.pi/180.*rlatlon[2])*np.sin(np.pi/180.*colatitude)
        xyz[1] = rlatlon[0]*np.sin(np.pi/180.*rlatlon[2])*np.sin(np.pi/180.*colatitude)
        xyz[2] = rlatlon[0]*np.cos(np.pi/180.*colatitude)
    else:
        raise ValueError('dimension of rlatlon should be 1 or 2')
    return xyz 

def polar2cart(rtheta):
    """
    Convert from polar to cartesian coordinates
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

def cart2polar(xy):
    """
    Convert from polar to cartesian coordinates
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
 
def getDestination(lat,lng,azimuth,distance):
    '''returns the lat an long of destination point 
    given the start lat, long, aziuth, and distance (in meters)'''
    R = constants.R #Radius of the Earth in m
    d = distance #Distance m 
    if lng > 180.: lng = lng -360.
    start = LatLon(lat,lng)
    end = start.destination(d,azimuth,R)
    return[end.lat, end.lon]

def calculateBearing(lat1,lng1,lat2,lng2):
    '''calculates the azimuth in degrees from start point to end point'''    
    if lng1 > 180.: lng1 = lng1 -360.
    start = LatLon(lat1,lng1)
    if lng2 > 180.: lng2 = lng2 -360.
    end = LatLon(lat2,lng2)
    bearing = start.initialBearingTo(end)
    return bearing

def getIntermediate(lat1,lng1,azimuth,distance,interval):
    '''returns every coordinate pair inbetween two coordinate 
    pairs given the desired interval. gcdelta and interval is great cirle dist in m'''
    R = constants.R #Radius of the Earth in m
    d = distance #Distance m 
    if lng1 > 180.: lng1 = lng1 -360.
    start = LatLon(lat1,lng1)
    end = start.destination(d,azimuth,R)
    steps = int(distance / interval)
    coords = []
    coords.append([lat1,lng1])
    for step in range(steps):
        counter = float(interval) * float(step+1)
        coord = getDestination(lat1,lng1,azimuth,counter)
        coords.append(coord)
    return coords
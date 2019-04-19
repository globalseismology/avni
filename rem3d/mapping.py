#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

from math import cos, pi, log, sin, tan, atan, atan2, sqrt, radians, degrees, asin, modf
import sys,os
import numpy as np #for numerical analysis
import multiprocessing
import codecs,json #printing output
from joblib import Parallel, delayed
import pdb    #for the debugger pdb.set_trace()
# from scipy.io import netcdf_file as netcdf #reading netcdf files
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
if (sys.version_info[:2] < (3, 0)): range = xrange

############################### PLOTTING ROUTINES ################################        
from .tools.trigd import atand,tand
from .tools.common import convert2nparray
from .f2py import ddelazgc # geolib library from NSW
from . import constants
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
        p1end = spher2cart(getDestinationLatLong(path1start[0],path1start[1],path1brngEnd,180.))
    else:
        path1def = 'endpoint' #path 1 defined by initial bearing
        p1end = path1brngEnd
    if isinstance(path2brngEnd, (float,int)):
        path2def = 'bearing'
        p2end = spher2cart(getDestinationLatLong(path2start[0],path2start[1],path2brngEnd,180.))
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
    return intersect,antipode

        
def midpoint(lat1, lon1, lat2, lon2):
    """Get the mid-point from positions in geographic coordinates.Input values as degrees"""

#Convert to radians
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)


    bx = cos(lat2) * cos(lon2 - lon1)
    by = cos(lat2) * sin(lon2 - lon1)
    lat3 = atan2(sin(lat1) + sin(lat2), \
           sqrt((cos(lat1) + bx) * (cos(lat1) \
           + bx) + by**2))
    lon3 = lon1 + atan2(by, cos(lat1) + bx)

    return [round(degrees(lat3), 2), round(degrees(lon3), 2)]


def get_distaz(eplat,eplon,stlat,stlon,num_cores=1):
    """Get the distance and azimuths from positions in geographic coordinates"""
    
    geoco = constants.geoco
    if isinstance(eplat,list): # if the input is a list loop 
        delta=[];azep=[];azst=[]
        # Standard checks on number of cores
        avail_cores = multiprocessing.cpu_count()
        if num_cores > avail_cores: 
            sys.exit("Number of cores requested ("+str(num_cores)+") is higher than available ("+str(avail_cores)+")")
        # Crate list of tuples of job arguments and pass to parallel using a helper routine     
        job_args = [(atand(geoco*tand(item_lat)),eplon[jj],atand(geoco*tand(stlat[jj])),stlon[jj]) for jj, item_lat in enumerate(eplat)]
        temp=Parallel(n_jobs=num_cores)(delayed(delazgc_helper)(ii) for ii in job_args)
        for il in temp: delta.append(il[0]);azep.append(il[1]);azst.append(il[2])

    elif isinstance(eplat,float):
        elat=atand(geoco*tand(eplat))
        elon=eplon
        slat=atand(geoco*tand(stlat))
        slon=stlon
        delta,azep,azst = ddelazgc(elat,elon,slat,slon)    
    else:    
        raise ValueError("get_distaz only takes list or floats")
    return delta,azep,azst

def delazgc_helper(args):
    return ddelazgc(*args)
    
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
 
def getDestinationLatLong(lat,lng,azimuth,distance):
    '''returns the lat an long of destination point 
    given the start lat, long, aziuth, and distance (in meters)'''
    R = constants.R #Radius of the Earth in m
    brng = radians(azimuth) #Bearing is degrees converted to radians.
    d = distance #Distance m 
    lat1 = radians(lat) #Current dd lat point converted to radians
    lon1 = radians(lng) #Current dd long point converted to radians
    lat2 = asin(sin(lat1) * cos(d/R) + cos(lat1)* sin(d/R)* cos(brng))
    lon2 = lon1 + atan2(sin(brng) * sin(d/R)* cos(lat1), cos(d/R)- sin(lat1)* sin(lat2))
    #convert back to degrees
    lat2 = degrees(lat2)
    lon2 = degrees(lon2)
    return[lat2, lon2]

def calculateBearing(lat1,lng1,lat2,lng2):
    '''calculates the azimuth in degrees from start point to end point'''
    startLat = radians(lat1)
    startLong = radians(lng1)
    endLat = radians(lat2)
    endLong = radians(lng2)
    dLong = endLong - startLong
    dPhi = log(tan(endLat/2.0+pi/4.0)/tan(startLat/2.0+pi/4.0))
    if abs(dLong) > pi:
         if dLong > 0.0:
             dLong = -(2.0 * pi - dLong)
         else:
             dLong = (2.0 * pi + dLong)
    bearing = (degrees(atan2(dLong, dPhi)) + 360.0) % 360.0;
    return bearing

def getintermediateLatLong(lat1,lng1,azimuth,gcdelta,interval):
    '''returns every coordinate pair inbetween two coordinate 
    pairs given the desired interval. gcdelta and interval is great cirle dist in m'''
#     d = getPathLength(lat1,lng1,lat2,lng2)
    remainder, dist = modf((gcdelta / interval))
    lat2,lng2=getDestinationLatLong(lat1,lng1,azimuth,gcdelta)
    counter = float(interval)
    coords = []
    coords.append([lat1,lng1])
    for distance in range(0,int(dist)):
        coord = getDestinationLatLong(lat1,lng1,azimuth,counter)
        counter = counter + float(interval)
        coords.append(coord)
    return coords

def interp_weights(xyz, uvw, d=3):
    """First, a call to sp.spatial.qhull.Dealunay is made to triangulate the irregular grid coordinates.
Then, for each point in the new grid, the triangulation is searched to find in which triangle (actually, in which simplex, which in your 3D case will be in which tetrahedron) does it lay.
The barycentric coordinates of each new grid point with respect to the vertices of the enclosing simplex are computed. From:
http://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
    """
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts, fill_value=np.nan):
    """An interpolated values is computed for that grid point, using the barycentric coordinates, and the values of the function at the vertices of the enclosing simplex. From:
    http://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids"""    
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value
    return ret
 

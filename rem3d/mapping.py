#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *


from math import cos, pi, log, sin, tan, atan, atan2, sqrt, radians, degrees, asin, modf
import sys,os
import numpy as np #for numerical analysis
# from scipy.io import netcdf_file as netcdf #reading netcdf files
import scipy.spatial.qhull as qhull

############################### PLOTTING ROUTINES ################################        
from . import constants
###############################

def atand(x):
    return degrees(atan(x))
    
def tand(x):
    return tan(radians(x))
        
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

    
def cart2spher(xyz):
    """Convert from cartesian to spherical coordinates
    http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node42.html
    """
    rlatlon = np.zeros(xyz.shape)
    xy = xyz[:,0]**2 + xyz[:,1]**2
    rlatlon[:,0] = np.sqrt(xy + xyz[:,2]**2)
#     ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
    rlatlon[:,1] = np.arctan2(xyz[:,2], np.sqrt(xy))/np.pi*180. # latitude
    rlatlon[:,2] = np.arctan2(xyz[:,1], xyz[:,0])/np.pi*180.
    return rlatlon
 
def spher2cart(rlatlon):
    """Convert from spherical to cartesian coordinates
    http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node42.html
    """
    xyz = np.zeros(rlatlon.shape)
    colatitude=90.-rlatlon[:,1]
    xyz[:,0] = rlatlon[:,0]*np.cos(np.pi/180.*rlatlon[:,2])*np.sin(np.pi/180.*colatitude)
    xyz[:,1] = rlatlon[:,0]*np.sin(np.pi/180.*rlatlon[:,2])*np.sin(np.pi/180.*colatitude)
    xyz[:,2] = rlatlon[:,0]*np.cos(np.pi/180.*colatitude)
    return xyz 

 
def getDestinationLatLong(lat,lng,azimuth,distance):
    '''returns the lat an long of destination point 
    given the start lat, long, aziuth, and distance'''
    R = 6378.1 #Radius of the Earth in km
    brng = radians(azimuth) #Bearing is degrees converted to radians.
    d = distance/1000 #Distance m converted to km
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

def getintermediateLatLong(lat1,lng1,azimuth,gcdeltakm,interval):
    '''returns every coordinate pair inbetween two coordinate 
    pairs given the desired interval. gcdeltakm and interval is great cirle dist in km'''
#     d = getPathLength(lat1,lng1,lat2,lng2)
    remainder, dist = modf((gcdeltakm / interval))
    lat2,lng2=getDestinationLatLong(lat1,lng1,azimuth,gcdeltakm)
    counter = float(interval)
    coords = []
    coords.append([lat1,lng1])
    for distance in xrange(0,int(dist)):
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


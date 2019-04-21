#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import numpy as np #for numerical analysis
# from scipy.io import netcdf_file as netcdf #reading netcdf files
import scipy.spatial.qhull as qhull

############################### PLOTTING ROUTINES ################################        
###############################
 
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
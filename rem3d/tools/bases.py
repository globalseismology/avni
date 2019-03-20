#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import pkgutil
import numpy as np
from collections import Counter

####################### IMPORT REM3D LIBRARIES  #######################################
from rem3d.f2py import vbspl,dbsplrem
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
        
    if isinstance(depths, (list,tuple,np.ndarray)):
        depths = np.asarray(depths)
    elif isinstance(depths, float):
        depths = np.asarray([depths])
    elif isinstance(depths, int):
        depths = np.asarray([float(depths)])
    else:
        raise TypeError('depths must be list or tuple, not %s' % type(depths))

    # find repeated values
    repeats = [item for item, count in Counter(knots).items() if count > 1]
    repeats_gt_2= [item for item, count in Counter(knots).items() if count > 2]
    if len(repeats_gt_2) != 0: raise ValueError('Cannot have more than 2 repetitions in knots')
    
    if len(repeats) > 0: # if there are repeated knots, splits it
        split_at = []
        for ii in range(len(repeats)):
            split_at.append(np.where(knots==repeats[ii])[0][1])
        knots_list = np.split(knots, split_at)
        for knots in knots_list: 
            if len(knots) < 4:
                raise ValueError('Atleast 4 knots need to be defined at or between '+str(min(knots))+' and '+str(max(knots))+' km') 
        
        jj = 0
        for depth in depths:
            jj = jj + 1
            for kk in range(len(knots_list)):
                # create the arrays as Fortran-contiguous
                splpts = np.array(knots_list[kk].tolist(), order = 'F')
                #Undefined if depth does not lie within the depth extents of knot points
                if depth < min(knots_list[kk]) or depth > max(knots_list[kk]): 
                    temp1 = temp2 = np.zeros_like(splpts)
                else:
                    (temp1, temp2) = vbspl(depth,len(splpts),splpts)
                if kk == 0:
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
        
    if isinstance(radius, (list,tuple,np.ndarray)):
        radius = np.asarray(radius)
    elif isinstance(radius, float):
        radius = np.asarray([radius])
    elif isinstance(radius, int):
        radius = np.asarray([float(radius)])
    else:
        raise TypeError('radius must be list or tuple, not %s' % type(depths))

    if len(radius_range) != 2 or not isinstance(radius_range, (list,tuple,np.ndarray)):
        raise TypeError('radius_range must be list , not %s' % type(radius_range))

    for irad in range(len(radius)):
        #Undefined if depth does not lie within the depth extents of knot points      
        if radius[irad] < min(radius_range) or radius[irad] > max(radius_range): 
            temp1 = temp2 = np.zeros(nsplines)
        else:
            (temp1, temp2) = dbsplrem(radius[irad],radius_range[0], radius_range[1],nsplines)
        if irad == 0:
            vercof = temp1; dvercof = temp2
        else:    
            vercof = np.vstack([vercof,temp1]) 
            dvercof = np.vstack([dvercof,temp2]) 
    return vercof, dvercof


def eval_polynomial(radius, radius_range, rnorm, types = ['CONSTANT','LINEAR']):
    """
    Evaluate the cubic spline know with second derivative as 0 at end points.
    
    Input parameters:
    ----------------
    
    radius: value or array of radii queried
    
    radius_range: limits of the radius limits of the region
    
    types: polynomial coefficients to be used for calculation. Options are : TOP,
                  TOP, BOTTOM, CONSTANT, LINEAR, QUADRATIC, CUBIC
    
    rnorm: normalization for radius, usually the radius of the planet
    
    Output:
    ------
    
    vercof : value of the polynomial coefficients at each depth, size (Nradius).
    
    """
        
    if isinstance(radius, (list,tuple,np.ndarray)):
        radiusin = np.asarray(radius)
    elif isinstance(radius, float):
        radiusin = np.asarray([radius])
    elif isinstance(radius, int):
        radiusin = np.asarray([float(radius)])
    else:
        raise TypeError('radius must be list or tuple, not %s' % type(radius))

    if len(radius_range) != 2 or not isinstance(radius_range, (list,tuple,np.ndarray)):
        raise TypeError('radius_range must be list , not %s' % type(radius_range))
    
    # keys in coefficients should be acceptable
    choices = ['TOP', 'BOTTOM', 'CONSTANT', 'LINEAR', 'QUADRATIC', 'CUBIC']
    assert(np.all([key in choices for key in types]))
    npoly = len(types)
    # first find whether CONSTANT/LINEAR or TOP/BOTTOM
    rbot=radius_range[0]/rnorm
    rtop=radius_range[1]/rnorm
    findtopbot = np.any([key in ['BOTTOM','TOP'] for key in types])
    findconstantlinear = np.any([key in ['CONSTANT','LINEAR'] for key in types])
    
    if findtopbot and findconstantlinear: raise ValueError('Cannot have both BOTTOM/TOP and CONSTANT/LINEAR as types in eval_polynomial')
    
    for irad in range(len(radiusin)):
        #Undefined if depth does not lie within the depth extents of knot points                      
        if radiusin[irad] < min(radius_range) or radiusin[irad] > max(radius_range): 
            temp = np.zeros(npoly)
            dtemp = np.zeros(npoly)
        else:
            rn=radiusin[irad]/rnorm  
            temp = np.zeros(npoly) 
            dtemp = np.zeros(npoly)
            for ii in range(npoly):
                if findconstantlinear:
                    if types[ii]=='CONSTANT':
                        temp[ii]=1.
                        dtemp[ii]=0.
                    elif types[ii]=='LINEAR':
                        temp[ii]=rn
                        dtemp[ii]=1.
                    elif types[ii]=='QUADRATIC':
                        temp[ii]=rn**2
                        dtemp[ii]=2.*rn
                    elif types[ii]=='CUBIC':
                        temp[ii]=rn**3
                        dtemp[ii]=3.*rn**2
                elif findtopbot:
                    if types[ii]=='TOP':
                        temp[ii] = 1.-(rn-rtop)/(rbot-rtop)
                        dtemp[ii]= -1./(rbot-rtop)
                    elif types[ii]=='BOTTOM':
                        temp[ii] = (rn-rtop)/(rbot-rtop)
                        dtemp[ii]= 1./(rbot-rtop)
                    elif types[ii]=='QUADRATIC':
                        temp[ii] = rn**2-rtop**2-(rn-rtop)*(rbot+rtop)
                        dtemp[ii]=2.*rn - 1./(rbot+rtop)
                    elif types[ii]=='CUBIC':
                        temp[ii]=rn**3-rtop**3-(rn-rtop)*(rbot**3-rtop**3)/(rbot-rtop)
                        dtemp[ii]=3.*rn**2 - 1.*(rbot**3-rtop**3)/(rbot-rtop)
        if irad == 0:
            vercof = temp; dvercof = dtemp
        else:    
            vercof = np.vstack([vercof,temp]) 
            dvercof = np.vstack([dvercof,dtemp]) 
    return vercof,dvercof

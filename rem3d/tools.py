#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import scipy.constants
import pkgutil
import os
import sys
import codecs,json #printing output
import numpy as np
from math import ceil
from collections import Counter
import pdb
import re
from configobj import ConfigObj

####################### IMPORT REM3D LIBRARIES  #######################################
from . import constants
from rem3d.f2py import vbspl,dbsplrem
#######################################################################################
    
def alphanum_key(s): 
    '''
    helper tool to sort lists in ascending numerical order (natural sorting),
    rather than lexicographic sorting
    '''
    return [int(c) if c.isdigit() else c for c in re.split('([0-9]+)', s)]


def krunge(n,x,h,y,f,m=0,phi=np.zeros(6),savey=np.zeros(6)):
    """
    some sort of integration or interpolation? x is 
    incremented on second and fourth calls. Resets itself after 
    5'th call.

    Input parameters:
    ----------------
    input: (these are guesses)
    m = call number
    n   = number of points
    f() = function evaluated at each point
    x   = independent variable
    h   = step size
    
    Output:
    ------
    y() = 
    """
    if len(y) > 6 or len(f) > 6: 
        raise ValueError ("len(y) > 6 or len(f) >  in krunge")
    m = m + 1    
    if m == 1:
        krunge=1
    elif m == 2:
        for j in np.arange(n): 
            savey[j] = y[j]
            phi[j]   = f[j]
            y[j] = savey[j]+(0.5*h*f[j])
        x = x + 0.5*h
        krunge = 1
    elif m == 3:
        for j in np.arange(n): 
            phi[j] = phi[j] + (2.0*f[j])
            y[j]   = savey[j] + (0.5*h*f[j])
        krunge = 1
    elif m == 4:
        for j in np.arange(n): 
            phi[j] = phi[j] + (2.0*f[j])
            y[j]   = savey[j] + (h*f[j])
        x = x + 0.5*h
        krunge = 1
    elif m == 5:
        for j in np.arange(n): 
            y[j] = savey[j] + (phi[j]+f[j])*h/6.0
        m    = 0
        krunge = 0
    return krunge,y,f,m,phi,savey
        

def firstnonspaceindex(string):
    """
    Gets the first non space index of a string
    """
    ifst=0
    ilst=len(string.rstrip('\n'))
    while string[ifst:ifst+1] == ' ' and ifst < ilst: ifst=ifst+1 
    if ilst-ifst <= 0: sys.exit("error reading model")
    return ifst,ilst

def get_fullpath(path):
    """
    Provides the full path by replacing . and ~ in path.
    """
    # Get the current directory    
    if path[0]=='.': path = os.path.dirname(os.path.abspath(__file__))+path[1:]   
    # If the path starts with tilde, replace with home directory
    if path[0]=='~': path=os.path.expanduser("~")+path[1:]
    return path
    
def listfolders(path):
    """
    Return a list of directories in a path
    """
    dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return dirs

def get_installdir(module='rem3d',checkwrite=True,checkenv=True):
    """
    Get the installation directory for any module. checkwrite checks for write access to the files.
    checkenv checks if the directory is specified as an environment variable.
    """
    installdir = None
    if checkenv:
        if os.environ.get(module+'_dir') is not None:
            installdir=os.environ.get(module+'_dir')
            print("Warning: Reading "+module+"_dir"+" from environment variables - "+installdir)

    if installdir is None:
        loader=pkgutil.find_loader(module)
        if loader is None:
            installdir = os.getcwd()
            print("Warning: installation directory not found for "+module+". Using current directory - "+installdir)
        else:
            installdir = os.path.dirname(loader.get_filename(module))
            if installdir is None:
                raise ValueError("Error: Specify "+module+"_dir environment variable (export "+module+"_dir=/path) or install "+module+" module.")
            if checkwrite:
                if not os.access(installdir, os.W_OK):
                    raise PermissionError("Error: Cannot I/O to "+module+" directory "+installdir+ " due to permissions. \
Specify "+module+"_dir environment variable or chdir to a different directory with I/O access.")
    return installdir

def get_filedir(module='rem3d',checkwrite=True,makedir=True):
    """
    Get the local files directory. Make a new directory if doesn't exist (makedir==True)
    """
    installdir = get_installdir(module=module,checkwrite=checkwrite)
    filedir = installdir+'/'+constants.localfilefolder
    if checkwrite and makedir: 
        if not os.path.exists(filedir):
            os.makedirs(filedir)        
    return filedir

def get_configdir(module='rem3d',checkwrite=True,makedir=True):
    """
    Get the directory containing configuration files. 
    Make a new directory if doesn't exist (makedir==True)
    """
    installdir = get_installdir(module=module,checkwrite=checkwrite)
    configdir = installdir+'/'+constants.configfolder
    if checkwrite and makedir: 
        if not os.path.exists(configdir):
            os.makedirs(configdir)        
    return configdir

def get_projections(checkwrite=True,makedir=True,type='radial'):
    """
    Get the file containing projection matrices. 
    Make a new directory if doesn't exist (makedir==True)
    """
    if type != 'radial' and type != 'lateral': 
        raise ValueError('type is undefined in get_projections')
    configdir = get_configdir(checkwrite=checkwrite,makedir=makedir)
    projections = configdir+'/projections.'+type+'.npz'
    exists = os.path.isfile(projections)
    return projections,exists

        
def writejson(nparray,filename,encoding='utf-8'):
    """Writes a json file from a numpy array"""
    
    listarray = nparray.tolist() # nested lists with same data, indices
    json.dump(listarray, codecs.open(filename, 'w', encoding=encoding), separators=(',', ':'), sort_keys=True, indent=4) ### this saves the array in .json format
    return

def readjson(filename,encoding='utf-8'):
    """Reading from a filename to a numpy array"""
        
    obj_text = codecs.open(filename, 'r', encoding=encoding).read()
    listarray = json.loads(obj_text)
    nparray = np.array(listarray)
    return nparray    

def uniquenumpyrow(a):
    """Gets the unique rows from a numpy array and the indices. e.g. to get unique lat-lon values"""
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    unique_a = a[idx]
    return unique_a,idx


def sanitised_input(prompt, type_=None, min_=None, max_=None, range_=None): 
    """Provide a user prompt with values between min-max or range of values e.g.
    For specific values:
    user_input = sanitised_input("Replace(r)/Ignore(i) this datum?", str.lower, range_=('r', 'i')
    For a range:
    age = sanitised_input("Enter your age: ", int, range_=xrange(100))
    """
    if min_ is not None and max_ is not None and max_ < min_: 
        raise ValueError("min_ must be less than or equal to max_.") 
    while True: 
        ui = input(prompt) 
        if type_ is not None: 
            try: 
                ui = type_(ui) 
            except ValueError: 
                print("Input type must be {0}.".format(type_.__name__)) 
                continue
        if max_ is not None and ui > max_: 
            print("Input must be less than or equal to {0}.".format(max_)) 
        elif min_ is not None and ui < min_: 
            print("Input must be greater than or equal to {0}.".format(min_)) 
        elif range_ is not None and ui not in range_: 
            if isinstance(range_, xrange): 
                template = "Input must be between {0} and {1}."
                print(template.format(range_[0],range_[-1])) 
            else: 
                template = "Input must be {0}."
                if len(range_) == 1: 
                    print(template.format(*range_)) 
                else: 
                    print(template.format(" or ".join((", ".join(map(str, 
                                                                     range_[:-1])), 
                                                       str(range_[-1]))))) 
        else: 
            return ui  
            
def getplanetconstants(planet = constants.planetpreferred, configfile = get_configdir()+'/'+constants.planetconstants):
    """
    Read the constants from configfile relevant to a planet to constants.py
    
    Input parameters:
    ----------------
    planet: planet option from configfile
    
    configfile: all the planet configurations are in this file. 
                Default option means read from tools.get_configdir()
    
    """
        
    if not os.path.isfile(configfile):
        raise IOError('No configuration file found: '+configfile)
    else:
        parser = ConfigObj(configfile)
    
    try:
        parser_select = parser[planet]
    except:
        raise IOError('No planet '+planet+' found in file '+configfile)
    constants.a_e = eval(parser_select['a_e']) # Equatorial radius
    constants.GM = eval(parser_select['GM']) # Geocentric gravitational constant m^3s^-2
    constants.G = eval(parser_select['G']) # Gravitational constant m^3kg^-1s^-2
    constants.f = eval(parser_select['f']) #flattening
    constants.omega = eval(parser_select['omega']) #Angular velocity in rad/s
    constants.M_true = eval(parser_select['M_true']) # Solid Earth mass in kg
    constants.I_true = eval(parser_select['I_true']) # Moment of inertia in m^2 kg
    constants.R = eval(parser_select['R']) # Radius of the Earth in m
    constants.rhobar = eval(parser_select['rhobar']) # Average density in kg/m^3
    # correction for geographic-geocentric conversion: 0.993277 for 1/f=297
    constants.geoco = (1.0 - constants.f)**2.  
   
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
    
    vercof, dvercof: value of the spline coefficients at each depth and its derivative.
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


def eval_polynomial(radius, radius_range, rnorm, coefficients = {'CONSTANT':0.,'LINEAR':0} ):
    """
    Evaluate the cubic spline know with second derivative as 0 at end points.
    
    Input parameters:
    ----------------
    
    radius: value or array of radii queried
    
    radius_range: limits of the radius limits of the region
    
    coefficients: polynomial coefficients to be used for calculation. Options are : TOP,
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
    
    # keys in coefficients should be acceptable
    choices = ['TOP', 'BOTTOM', 'CONSTANT', 'LINEAR', 'QUADRATIC', 'CUBIC']
    assert(np.all([key in choices for key in coefficients.keys()]))
    npoly = len(coefficients.keys())
    rfnval = {}
    for choice in choices:
        if choice in coefficients: #if the coefficient is defined, store in rfnval
            rfnval[choice] = coefficients[choice]
        else:
            rfnval[choice] = 0.
    # firstfine CONSTANT and linear from TOP and BOTTOM
    rbot=radius_range[0]/rnorm
    rtop=radius_range[1]/rnorm
    findtopbot = np.any([key in ['BOTTOM','TOP'] for key in coefficients.keys()])
    if findtopbot:
        s1=rfnval['BOTTOM']-rfnval['TOP']-rfnval['QUADRATIC']*(rbot**2-rtop**2)-rfnval['CUBIC']*(rbot**3-rtop**3)
        s1=s1/(rbot-rtop)
        rfnval['LINEAR']=s1
        rfnval['CONSTANT']=rfnval['TOP']-rfnval['LINEAR']*rtop- rfnval['QUADRATIC']*rtop**2-rfnval['CUBIC']*rtop**3

    if len(radius_range) != 2 or not isinstance(radius_range, (list,tuple,np.ndarray)):
        raise TypeError('radius_range must be list , not %s' % type(radius_range))

    for irad in range(len(radiusin)):
        #Undefined if depth does not lie within the depth extents of knot points                      
        if radiusin[irad] < min(radius_range) or radiusin[irad] > max(radius_range): 
            temp = 0.
        else:
            rn=radiusin[irad]/rnorm                
            temp=rfnval['CONSTANT']+rfnval['LINEAR']*rn+rfnval['QUADRATIC']*rn**2+ rfnval['CUBIC']*rn**3
        
        if irad == 0:
            if len(radiusin) == 1:
                vercof = temp
            else:
                vercof = [temp]
        else:    
            vercof.append(temp)
    # convert to appropriate type based on input
    if isinstance(radius, np.ndarray):
        vercof = np.array(vercof)
    elif isinstance(radius, tuple):
        vercof = tuple(vercof)
    return vercof

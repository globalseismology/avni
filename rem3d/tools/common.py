#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int,list,tuple

import pkgutil
import os
import codecs,json #printing output
import numpy as np
import re
from configobj import ConfigObj
from six import string_types # to check if variable is string using isinstance
import ntpath
import ast

####################### IMPORT REM3D LIBRARIES  #######################################
from .. import constants
#######################################################################################

def stage(file,overwrite=False):
    """
    Stages a file in the rem3d file directories for testing
    """
    filedir = get_filedir()
    if not os.path.isfile(file): raise IOError(file+' not found')
    outlink = filedir+'/'+ntpath.basename(file)
    try:
        os.symlink(file, outlink)
    except OSError:
        if overwrite:
            os.unlink(outlink)
            os.symlink(file, outlink)
        else:
            print('Warning: a link to file '+ntpath.basename(file)+' exists within REM3D. Use overwrite=True to overwrite the staged link.')
    return

def convert2nparray(value,int2float = True):
    """
    Converts input value to a float numpy array

    int2float: convert integer to floats, if true
    """
    if isinstance(value, (list,tuple,np.ndarray)):
        outvalue = np.asarray(value)
    elif isinstance(value, float):
        outvalue = np.asarray([value])
    elif isinstance(value, (int,np.int64)):
        if int2float:
            outvalue = np.asarray([float(value)])
        else:
            outvalue = np.asarray([value])
    elif isinstance(value,string_types):
        outvalue = np.asarray([value])
    else:
        raise TypeError('input must be list or tuple, not %s' % type(value))
    return outvalue


def precision_and_scale(x):
    """
    Returns precision and scale of a float
    """
    max_digits = 14
    int_part = np.int(abs(x))
    magnitude = 1 if int_part == 0 else np.int(np.log10(int_part)) + 1
    if magnitude >= max_digits:
        return (magnitude, 0)
    frac_part = abs(x) - int_part
    multiplier = 10 ** (max_digits - magnitude)
    frac_digits = multiplier + np.int(multiplier * frac_part + 0.5)
    while frac_digits % 10 == 0:
        frac_digits /= 10
    scale = np.int(np.log10(frac_digits))
    return (magnitude + scale, scale)

def alphanum_key(s):
    '''
    helper tool to sort lists in ascending numerical order (natural sorting),
    rather than lexicographic sorting
    '''
    return [int(c) if c.isdigit() else c for c in re.split('([0-9]+)', s)]

def diffdict(first_dict,second_dict):
    '''
    helper tool to get difference in two dictionaries
    '''
    return { k : second_dict[k] for k in set(second_dict) - set(first_dict) }

def df2nparray(dataframe):
    '''
    helper tool to return the named numpy array of the pandas dataframe
    '''
    columns = {}
    for ii in range(dataframe.shape[1]): columns[ii] = str(ii)
    dataframe.rename(columns = columns,inplace=True)
    ra = dataframe.to_records(index=False)
    return np.asarray(ra)

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
    if ilst-ifst <= 0: raise ValueError("error reading model")
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

def get_cptdir(module='rem3d',checkwrite=True,makedir=True):
    """
    Get the directory with color palettes. Make a new directory if doesn't exist (makedir==True)
    """
    filedir = get_filedir(module=module,checkwrite=checkwrite,makedir=makedir)
    cptdir = filedir+'/'+constants.cptfolder
    if checkwrite and makedir:
        if not os.path.exists(cptdir):
            os.makedirs(cptdir)
    return cptdir

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
    constants.a_e = ast.literal_eval(parser_select['a_e']) # Equatorial radius
    constants.GM = ast.literal_eval(parser_select['GM']) # Geocentric gravitational constant m^3s^-2
    constants.G = ast.literal_eval(parser_select['G']) # Gravitational constant m^3kg^-1s^-2
    try:
        constants.f = ast.literal_eval(parser_select['f']) #flattening
    except KeyError:
        try:
            constants.f = 1./ast.literal_eval(parser_select['1/f']) #flattening
        except:
            raise KeyError('need either flattening (f) or inverse flattening (1/f) for '+planet+' in '+configfile)
    constants.omega = ast.literal_eval(parser_select['omega']) #Angular velocity in rad/s
    constants.M_true = ast.literal_eval(parser_select['M_true']) # Solid Earth mass in kg
    constants.I_true = ast.literal_eval(parser_select['I_true']) # Moment of inertia in m^2 kg
    constants.R = ast.literal_eval(parser_select['R']) # Radius of the Earth in m
    constants.rhobar = ast.literal_eval(parser_select['rhobar']) # Average density in kg/m^3
    constants.deg2km = ast.literal_eval(parser_select['deg2km']) #length of 1 degree in km
    constants.deg2m = constants.deg2km * 1000. #length of 1 degree in m
    # correction for geographic-geocentric conversion: 0.993277 for 1/f=297
    try:
        print('... Re - Initialized rem3d module with constants for '+planet+' from '+parser_select['cite']+' from geocentric correction '+str(constants.geoco))
        constants.geoco = (1.0 - constants.f)**2.
    except AttributeError:
        constants.geoco = (1.0 - constants.f)**2.


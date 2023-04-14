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
import pint # For SI units
import decimal
from numba import jit
import pdb

####################### IMPORT AVNI LIBRARIES  #######################################
from .. import constants
#######################################################################################

def stage(file,overwrite=False):
    """
    Stages a file in the avni file directories for testing
    """
    filedir = get_filedir() #AVNI file directory
    stagedfile = get_fullpath(file)
    if not os.path.isfile(stagedfile): raise IOError(stagedfile+' not found')
    outlink = filedir+'/'+ntpath.basename(file)
    try:
        os.symlink(stagedfile, outlink)
    except OSError:
        if overwrite and os.path.islink(outlink):
            os.unlink(outlink)
            os.symlink(stagedfile, outlink)
            print('WARNING: overwriting an existing staged link named '+ntpath.basename(file))
        else:
            raise IOError('A link to an actual file '+ntpath.basename(file)+' (not symlink) exists within AVNI. Delete the file '+outlink+' first before proceeding.')
    return

def parse_line(line,rx_dict):
    """
    Function used to parse line with key word from rx_dict

    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex
    """
    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None

@jit(nopython=True)
def ifwithindepth(start_depths,end_depths,depth_in_km):
    if depth_in_km.ndim != 1: raise ValueError('only 1-D array depth_in_km allowed')
    output = np.ones_like(depth_in_km,dtype=np.int64)
    output[:] = -1
    for ii,depth in enumerate(depth_in_km):
        for jj,start in enumerate(start_depths):
            end = end_depths[jj]
            if depth >= start and depth <= end:
                output[ii]=jj
                break
    return output

def makegrid(latitude,longitude,depth_in_km=None):
    """
    Make a 2D or 3D grid out of input locations and depths.

    grid: make a grid by unraveling (depth_in_km,latitude,longitude)
    """

    # convert to numpy arrays
    latitude = convert2nparray(latitude)
    longitude = convert2nparray(longitude)
    nlat = len(latitude)
    nlon = len(longitude)

    #checks
    if depth_in_km==None:
        if not (latitude.ndim == longitude.ndim == 1): raise ValueError("latitude, longitude or depth_in_km should be one-dimensional arrays")
        nrows = nlat*nlon
        longitude,latitude = np.meshgrid(longitude,latitude)
        longitude = longitude.ravel()
        latitude = latitude.ravel()
        if not (len(latitude) == len(longitude)): raise ValueError("latitude, longitude or depth_in_km should be of same length if not making grid = False")
        return nrows,latitude,longitude
    else:
        if not (latitude.ndim == longitude.ndim == depth_in_km.ndim == 1): raise ValueError("latitude, longitude or depth_in_km should be one-dimensional arrays")
        depth_in_km = convert2nparray(depth_in_km)
        ndep = len(depth_in_km)
        nrows = ndep*nlat*nlon
        depth_tmp = np.zeros(nrows)
        for indx,depth in enumerate(depth_in_km):
            depth_tmp[indx*nlat*nlon:(indx+1)*nlat*nlon] = depth * np.ones(nlat*nlon)
            latitude = np.tile(latitude,len(depth_in_km))
            longitude = np.tile(longitude,len(depth_in_km))
            depth_in_km = depth_tmp
        if not (len(latitude) == len(longitude) == len(depth_in_km)): raise ValueError("latitude, longitude or depth_in_km should be of same length if not making grid = False")
        return nrows,latitude,longitude,depth_in_km

def convert2nparray(value,int2float = True,allowstrings=True):
    """
    Converts input value to a float numpy array. Boolean are returned as Boolean arrays.

    int2float: convert integer to floats, if true

    stringallowed: check if value has strings
    """
    if isinstance(value, (list,tuple,np.ndarray)):
        outvalue = np.asarray(value)
    elif isinstance(value, bool):
        outvalue = np.asarray([value])
    elif isinstance(value, float):
        outvalue = np.asarray([value])
    elif isinstance(value, (int,np.int64)):
        if int2float:
            outvalue = np.asarray([float(value)])
        else:
            outvalue = np.asarray([value])
    elif isinstance(value,string_types):
        if not allowstrings: raise TypeError('input cannot be a string')
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
    '''helper tool to get difference in two dictionaries
    '''
    return { k : second_dict[k] for k in set(second_dict) - set(first_dict) }

def equaldict(first_dict,second_dict):
    '''helper tool to check if two dictionaries are equal
    '''
    #checks=[]
    #for k in set(realization.metadata):
    #    checks.extend(convert2nparray(first_dict==second_dict))
    #return np.all(checks)

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
    if path[0]=='.': path = os.path.abspath(path)
    if len(path) > 2 :
        if path[:2]=='..': path = os.path.abspath(path)
    # If the path starts with tilde, replace with home directory
    if path[0]=='~': path = os.path.expanduser(path)
    # if the expanded path does not have / as first character it is from current directory
    # if only file name provided append current directory
    if ntpath.basename(path) == path or path[0] != '/': path = os.path.abspath('./'+path)
    return path

def listfolders(path):
    """
    Return a list of directories in a path
    """
    dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return dirs

def get_installdir(module='avni',checkwrite=True,checkenv=True):
    """
    Get the installation directory for any module. checkwrite checks for write access to the files.
    checkenv checks if the directory is specified as an environment variable.
    """
    installdir = None
    if checkenv:
        if os.environ.get(module+'_dir') is not None:
            installdir=os.environ.get(module+'_dir')
            # print("Warning: Reading "+module+"_dir"+" from environment variables - "+installdir)

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

def get_filedir(module='avni',subdirectory=None,checkwrite=True,makedir=True):
    """
    Get the local files directory. Make a new directory if doesn't exist (makedir==True)
    """
    installdir = get_installdir(module=module,checkwrite=checkwrite)
    filedir = os.path.join(installdir,constants.localfilefolder)
    if subdirectory is not None: filedir = os.path.join(filedir,subdirectory)

    if checkwrite and makedir:
        if not os.path.exists(filedir): os.makedirs(filedir)
    return filedir

def get_cptdir(module='avni',checkwrite=True,makedir=True):
    """
    Get the directory with color palettes. Make a new directory if doesn't exist (makedir==True)
    """
    return get_filedir(module=module,subdirectory=constants.cptfolder,
                        checkwrite=checkwrite,makedir=makedir)

def get_configdir(module='avni',checkwrite=True,makedir=True):
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

def get_projections(checkwrite=True,makedir=True,types='radial'):
    """
    Get the file containing projection matrices.
    Make a new directory if doesn't exist (makedir==True)
    """
    if types != 'radial' and types != 'lateral':
        raise ValueError('types is undefined in get_projections')
    configdir = get_configdir(checkwrite=checkwrite,makedir=makedir)
    projections = configdir+'/projections.'+types+'.npz'
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

def appendunits(ureg=constants.ureg,system='mks',unitsfile = get_configdir()+'/'+constants.customunits):
    """
    Append the custom units from unitsfile to ureg registry

    Input parameters:
    ----------------
    ureg: input unit registry. If None, initialize within to system

    system: default unit system. If not the same as ureg, change it.

    unitsfile: additional definitions to add to ureg
    """
    if ureg is None:
        ureg = pint.UnitRegistry(system=system)
    else:
        if ureg.default_system != system: ureg.default_system = system
    ureg.load_definitions(unitsfile)
    constants.ureg = ureg

def convert2units(valstring):
    """
    Returns the value with units. Only space allowed is that between value and unit.
    """
    vals = valstring.split()
    if len(vals) == 1: #if no unit is provided
        return eval(vals[0])*constants.ureg('dimensionless')
    elif len(vals) == 2: # first is value, second unit
        return eval(vals[0])*constants.ureg(vals[1])
    else:
        raise ValueError('only space allowed is that between value and unit')

def decimals(value):
    """
    Returns the number of decimals in a float
    """
    value = convert2nparray(value)
    output = np.zeros(len(value), dtype=int)
    for indx, val in enumerate(value):
        d = decimal.Decimal(str(val))
        output[indx] = abs(d.as_tuple().exponent)
    return output if len(output) > 1 else output[0]
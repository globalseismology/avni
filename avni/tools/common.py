#!/usr/bin/env python

#####################  IMPORT STANDARD MODULES   #########################

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
from six import string_types # to check if variable is string using isinstance
import ntpath
import pint # For SI units
import decimal
from numba import jit
import typing as tp
import pandas as pd
import warnings

####################### IMPORT AVNI LIBRARIES  ###########################

from .. import constants

##########################################################################

def stage(file: str, overwrite: bool = False):
    """Stages a file in the AVNI file directories for testing

    Parameters
    ----------
    file : str
        file name to be staged
    overwrite : bool, optional
        overwite any existing symlink or file, by default False

    Raises
    ------
    IOError
        File not found or a symlink already exists

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
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
            warnings.warn('WARNING: overwriting an existing staged link named '+ntpath.basename(file))
        else:
            raise IOError('A link to an actual file '+ntpath.basename(file)+' (not symlink) exists within AVNI. Delete the file '+outlink+' first before proceeding.')
    return

def parse_line(line: str,rx_dict: dict):
    """Function used to parse line with key word from rx_dict

    Parameters
    ----------
    line : str
        Line to search
    rx_dict : dict
        Do a regex search against all defined regexes

    Returns
    -------
    key,match
        Result of the first matching regex, None if not found

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None

@jit(nopython=True)
def ifwithindepth(start_depths: np.ndarray,end_depths: np.ndarray,depth_in_km: np.ndarray):
    """Check if a set of depths are within a range of depths

    Parameters
    ----------
    start_depths : np.ndarray
        Starting depth for the range
    end_depths : np.ndarray
        End depth for the range
    depth_in_km : np.ndarray
        Depths to check for

    Returns
    -------
    output
        index of depth range that each `depth_in_km` belongs to

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
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

def makegrid(latitude: tp.Union[list,tuple,np.ndarray],
             longitude: tp.Union[list,tuple,np.ndarray],
             depth_in_km: tp.Union[None,list,tuple,np.ndarray] = None):
    """Make a 2D or 3D grid out of input locations and depths.

    Parameters
    ----------
    latitude : tp.Union[list,tuple,np.ndarray]
        Latitudes of locations queried
    longitude : tp.Union[list,tuple,np.ndarray]
        Longitudes of locations queried
    depth_in_km : tp.Union[None,list,tuple,np.ndarray], optional
        Depths of locations queried, by default None

    Returns
    -------
    nrows,latitude,longitude,[optional: depth_in_km]
        A 2D or 3D grid with `nrows` rows found by unraveling (depth_in_km,latitude,longitude)

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
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

def convert2nparray(value: tp.Union[str,list,tuple,float,np.int64,np.ndarray,bool],
                    int2float: bool = True,
                    allowstrings: bool = True) -> np.ndarray:
    """Converts input value to a float numpy array. Boolean are returned as Boolean arrays.

    Parameters
    ----------
    value : tp.Union[str,list,tuple,float,np.int64,np.ndarray,bool]
        A single value or a set of values.
    int2float : bool, optional
        Convert integer to floats, by default True
    allowstrings : bool, optional
        Check if value has strings, by default True

    Returns
    -------
    np.ndarray
        Output values as a numpy array.

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # checks
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


def precision_and_scale(x: float, max_digits: int = 14) -> tp.Tuple[int, int]:
    """Returns precision and scale of a float

    Parameters
    ----------
    x : float
        A floating point number
    max_digits : int, optional
        Maximum digits or magnitude, by default 14

    Returns
    -------
    tp.Tuple[int, int]
        Returns precision and scale of a float

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

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

def alphanum_key(s: str) -> list:
    """Helper tool to sort lists in ascending numerical order (natural sorting),
    rather than lexicographic sorting

    Parameters
    ----------
    s : str
        string to sort

    Returns
    -------
    list
        Sorted list

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    return [int(c) if c.isdigit() else c for c in re.split('([0-9]+)', s)]

def diffdict(first_dict: dict,second_dict: dict) -> dict:
    """Helper tool to get difference in two dictionaries

    Parameters
    ----------
    first_dict : dict
        First dictionary
    second_dict : dict
        Second dictionary

    Returns
    -------
    dict
        Difference

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    return { k : second_dict[k] for k in set(second_dict) - set(first_dict) }

#def equaldict(first_dict,second_dict):
#    '''helper tool to check if two dictionaries are equal
#    '''
    #checks=[]
    #for k in set(realization.metadata):
    #    checks.extend(convert2nparray(first_dict==second_dict))
    #return np.all(checks)

def df2nparray(dataframe: pd.DataFrame) -> np.ndarray:
    """Helper tool to return the named numpy array of the pandas dataframe

    Parameters
    ----------
    dataframe : pd.DataFrame
        Input Pandas Dataframe

    Returns
    -------
    np.ndarray
        Equivalent named numpy array

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    columns = {}
    for ii in range(dataframe.shape[1]): columns[ii] = str(ii)
    dataframe.rename(columns = columns,inplace=True)
    ra = dataframe.to_records(index=False)
    return np.asarray(ra)

def krunge(n:int,x: float,
           h: tp.Union[int,float],
           y: tp.Union[list,tuple,np.ndarray],
           f: tp.Union[list,tuple,np.ndarray],
           m: int = 0,
           phi: np.ndarray = np.zeros(6),
           savey: np.ndarray = np.zeros(6)):
    """Some sort of integration or interpolation?
    x is incremented on second and fourth calls. Resets itself after 5'th call.

    Parameters
    ----------
    n : int
        number of points
    x : float
        indepent variable for points
    h : tp.Union[int,float]
        step size
    y : tp.Union[list,tuple,np.ndarray]
        function evaluated at each point
    f : tp.Union[list,tuple,np.ndarray]
        function evaluated at each point
    m : int, optional
        call number, by default 0
    phi : np.ndarray, optional
        Intermediate variable, by default np.zeros(6)
    savey : np.ndarray, optional
        Intermediate variable, by default np.zeros(6)
    """

    raise NotImplementedError('This function has not been tested against Fortran codes')

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


def firstnonspaceindex(string: str) -> tp.Tuple[int,int]:
    """Gets the first and last non-space index of a string

    Parameters
    ----------
    string : str
        String to search

    Returns
    -------
    tp.Tuple[int,int]
        First and last non-space index

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    ifst=0
    ilst=len(string.rstrip('\n'))
    while string[ifst:ifst+1] == ' ' and ifst < ilst: ifst=ifst+1
    if ilst-ifst <= 0: raise ValueError("error reading model")
    return ifst,ilst

def get_fullpath(path: str) -> str:
    """Provides the full path by replacing . and ~ in path

    Parameters
    ----------
    path : str
        Input file or folder path

    Returns
    -------
    str
        Full path after replacing . and ~ in path

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
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

def listfolders(path: str) -> list:
    """Return a list of directories in a path

    Parameters
    ----------
    path : str
        Input file or folder path

    Returns
    -------
    list
        List of directories

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return dirs

def get_installdir(module: str = 'avni',
                   checkwrite: bool = True,
                   checkenv: bool =True) -> str:
    """Get the installation directory for any module in the current machine.

    Parameters
    ----------
    module : str, optional
        Module to search for, by default 'avni'
    checkwrite : bool, optional
        Checks for write access to the files, by default True
    checkenv : bool, optional
        Checks if the directory is specified as an environment variable, by default True

    Returns
    -------
    str
        Path to installation directory

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
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

def get_filedir(module: str = 'avni',
                subdirectory: tp.Union[None,str] = None,
                checkwrite: bool = True,
                makedir: bool = True) -> str:
    """Get the local files directory for any module in the current machine.

    Parameters
    ----------
    module : str, optional
        Module to search for, by default 'avni'
    subdirectory : tp.Union[None,str], optional
        A directory inside the main directory, by default None
    checkwrite : bool, optional
        Checks for write access to the files, by default True
    makedir : bool, optional
        Make a new directory if it doesn't exist, by default True

    Returns
    -------
    str
        Path to local files directory

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    installdir = get_installdir(module=module,checkwrite=checkwrite)
    filedir = os.path.join(installdir,constants.localfilefolder)
    if subdirectory is not None: filedir = os.path.join(filedir,subdirectory)

    if checkwrite and makedir:
        if not os.path.exists(filedir): os.makedirs(filedir)
    return filedir

def get_cptdir(module: str = 'avni',
               checkwrite: bool = True,
               makedir: bool = True) -> str:
    """Get the directory with color palettes for any module in the current machine.

    Parameters
    ----------
    module : str, optional
        Module to search for, by default 'avni'
    checkwrite : bool, optional
        Checks for write access to the files, by default True
    makedir : bool, optional
        Make a new directory if it doesn't exist, by default True

    Returns
    -------
    str
        Path to local CPT color palette directory

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    return get_filedir(module=module,subdirectory=constants.cptfolder,
                        checkwrite=checkwrite,makedir=makedir)

def get_configdir(module: str = 'avni',
                  checkwrite: bool = True,
                  makedir: bool = True) -> str:
    """Get the directory containing configuration files.

    Parameters
    ----------
    module : str, optional
        Module to search for, by default 'avni'
    checkwrite : bool, optional
        Checks for write access to the files, by default True
    makedir : bool, optional
        Make a new directory if it doesn't exist, by default True

    Returns
    -------
    str
        Path to local config directory as specified in :py:func:`constants.configfolder`

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    installdir = get_installdir(module=module,checkwrite=checkwrite)
    configdir = installdir+'/'+constants.configfolder
    if checkwrite and makedir:
        if not os.path.exists(configdir):
            os.makedirs(configdir)
    return configdir

def get_projections(checkwrite: bool = True,
                    makedir: bool = True,
                    types: str = 'radial') -> tp.Tuple[str,bool]:
    """Get the file containing projection matrices.

    Parameters
    ----------
    checkwrite : bool, optional
        Checks for write access to the files, by default True
    makedir : bool, optional
        Make a new directory if doesn't exist, by default True
    types : str, optional
        Type can be radial or lateral, by default 'radial'

    Returns
    -------
    tp.Tuple[str,bool]
        First entry is location of projection file
        Second end is whether this file already exists in the file system

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    if types != 'radial' and types != 'lateral':
        raise ValueError('types is undefined in get_projections')
    configdir = get_configdir(checkwrite=checkwrite,makedir=makedir)
    projections = configdir+'/projections.'+types+'.npz'
    exists = os.path.isfile(projections)
    return projections,exists

def writejson(nparray: np.ndarray,
              filename: str,
              encoding: str = 'utf-8'):
    """Writes a json file from a numpy array

    Parameters
    ----------
    nparray : np.ndarray
        Input numpy array to write to JSON
    filename : str
        Output JSON file name
    encoding : str, optional
        File encoding, by default 'utf-8'

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    listarray = nparray.tolist() # nested lists with same data, indices
    json.dump(listarray, codecs.open(filename, 'w', encoding=encoding), separators=(',', ':'), sort_keys=True, indent=4) ### this saves the array in .json format
    return

def readjson(filename: str,encoding: str = 'utf-8') -> np.ndarray:
    """Reading from a JSON file to a numpy array

    Parameters
    ----------
    filename : str
        Input JSON file name
    encoding : str, optional
        File encoding, by default 'utf-8'

    Returns
    -------
    np.ndarray
        Output numpy array

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    obj_text = codecs.open(filename, 'r', encoding=encoding).read()
    listarray = json.loads(obj_text)
    nparray = np.array(listarray)
    return nparray

def uniquenumpyrow(a: np.ndarray):
    """Gets the unique rows from a numpy array and the indices. e.g. to get unique lat-lon values

    Parameters
    ----------
    a : np.ndarray
        A numpy aray with multiple rows that needs to searched

    Returns
    -------
    unique_a,idx
        Unique rows and their indices

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    unique_a = a[idx]
    return unique_a,idx


def sanitised_input(prompt: str, type_=None, min_=None, max_=None, range_=None):
    """Provide a user prompt with values between min-max or range of values

    For specific values:
    user_input = sanitised_input("Replace(r)/Ignore(i) this datum?", str.lower, range_=('r', 'i')
    For a range:
    age = sanitised_input("Enter your age: ", int, range_=xrange(100))

    Parameters
    ----------
    prompt : str
        Prompt to the user
    type_ : optional
        type of input needed, by default None
    min_ : optional
        Minimum value, by default None
    max_ : optional
        Maximum value, by default None
    range_ : optional
        Range of values, by default None

    Returns
    -------
    ui
        User input

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
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
            if isinstance(range_, range):
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

def appendunits(ureg=constants.ureg,
                system: str = 'mks',
                unitsfile: str = get_configdir()+'/'+constants.customunits):
    """Append the custom units from unitsfile to `ureg` registry

    Parameters
    ----------
    ureg : _type_, optional
        Input unit registry, by default constants.ureg
    system : str, optional
        Default unit system; if not the same as ureg changes it, by default 'mks'
    unitsfile : str, optional
        additional definitions to add to `ureg`, by default get_configdir()+'/'+constants.customunits

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    if ureg is None:
        ureg = pint.UnitRegistry(system=system)
    else:
        if ureg.default_system != system: ureg.default_system = system
    ureg.load_definitions(unitsfile)
    constants.ureg = ureg

def convert2units(valstring: str):
    """Returns the value with units from a string that has units.

    Parameters
    ----------
    valstring : str
        A string to convert. Only space allowed is that between value and unit.

    Returns
    -------
    vals
        Value with units

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    vals = valstring.split()
    if len(vals) == 1: #if no unit is provided
        return eval(vals[0])*constants.ureg('dimensionless')
    elif len(vals) == 2: # first is value, second unit
        return eval(vals[0])*constants.ureg(vals[1])
    else:
        raise ValueError('only space allowed is that between value and unit')

def decimals(value: tp.Union[list,tuple,np.ndarray]) -> np.ndarray:
    """Returns the number of decimals in a set of floating point numbers.

    Parameters
    ----------
    value : tp.Union[list,tuple,np.ndarray]
        Set of floating point numbers

    Returns
    -------
    np.ndarray
        Number of decimals in each float

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    value = convert2nparray(value)
    output = np.zeros(len(value), dtype=int)
    for indx, val in enumerate(value):
        d = decimal.Decimal(str(val))
        output[indx] = abs(d.as_tuple().exponent)
    return output if len(output) > 1 else output[0]
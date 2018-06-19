from __future__ import absolute_import
import scipy.constants
import pkgutil
import os
import sys
import codecs,json #printing output
import numpy as np
from math import ceil
####################### IMPORT REM3D LIBRARIES  #######################################
from . import constants
#######################################################################################

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

def get_installdir(checkwrite=True,checkenv=True):
    """
    Get the installation directory for the rem3d module. checkwrite checks for write access to the files.
    checkenv checks if the directory is specified as an environment variable.
    """
    installdir = None
    if checkenv:
        if os.environ.get('REM3Ddir') is not None:
            installdir=os.environ.get('REM3Ddir')
            print("Warning: Reading REM3Ddir from environment variables - "+installdir)

    if installdir is None:
        loader=pkgutil.find_loader('rem3d')
        if loader is None:
            installdir = os.getcwd()
            print("Warning: installation directory not found for rem3d. Using current directory - "+installdir)
        else:
            installdir = os.path.dirname(loader.get_filename('rem3d'))
            if installdir is None:
                raise ValueError("Error: Specify REM3Ddir environment variable (export REM3Ddir=/path/to/rem3d/) or install rem3d module. ")
            if checkwrite:
                if not os.access(installdir, os.W_OK):
                    raise PermissionError("Error: Cannot I/O to rem3d directory "+installdir+ " due to permissions. \
Specify rem3d_dir environment variable or chdir to a different directory with I/O access. ")
    return installdir

def get_filedir(checkwrite=True,makedir=True):
    """
    Get the local files directory. Make a new directory if doesn't exist (makedir==True)
    """
    installdir = get_installdir(checkwrite=checkwrite)
    filedir = installdir+'/'+constants.localfilefolder
    if checkwrite and makedir: 
        if not os.path.exists(filedir):
            os.makedirs(filedir)        
    return filedir

def get_configdir(checkwrite=True,makedir=True):
    """
    Get the directory containing configuration files. 
    Make a new directory if doesn't exist (makedir==True)
    """
    installdir = get_installdir(checkwrite=checkwrite)
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

from __future__ import absolute_import
import scipy.constants
import ConfigParser
import pdb    #for the debugger pdb.set_trace()
import pkgutil
import os
import sys
import codecs,json #printing output
import numpy as np
####################### IMPORT REM3D LIBRARIES  #######################################
from . import constants
#######################################################################################

def get_fullpath(path):
    """Provides the full path by replacing . and ~ in path."""
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
    installdir='not_set'
    ierror=0
    if checkenv: 
        if os.environ.get('REM3Ddir') is not None: 
            installdir=os.environ.get('REM3Ddir')
            print "Warning: Reading REM3Ddir from environment variables. "
    
    if installdir == 'not_set':        
        loader=pkgutil.find_loader('rem3d')
        if loader is None:
            print "Warning: installation directory not found for rem3d. Using current directory "
            installdir = os.getcwd()
        else:
            installdir = os.path.dirname(loader.get_filename('rem3d'))
            if installdir is None: 
                print "Error: Specify rem3d_dir environment variable or install rem3d module. "
                ierror=1
            if checkwrite:
                if not os.access(installdir, os.W_OK):
                    print "Error: Cannot I/O to rem3d directory "+installdir+ " due to permissions. \
Specify rem3d_dir environment variable or chdir to a different directory with I/O access. "
                    ierror=2
    return installdir, ierror
    
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


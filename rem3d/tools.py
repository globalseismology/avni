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
        ui = raw_input(prompt) 
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

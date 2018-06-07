#!/usr/bin/env python
"""
This script/module contains routines that are used to analyze Earth models and files that
contain them.
"""

#####################  IMPORT STANDARD MODULES   ######################################   

import os,sys
import numpy as np #for numerical analysis

#####################

def readepixfile(filename):
    """Read .epix file format from a file. 
    
    Parameters
    ----------
    
    filename : Name of the file containing four columns
              (latitude, longitude, pixel_size, value)
    
    """

    currentdir=os.getcwd()
    try: 
        f = open(filename, 'r')
        epixarr=np.genfromtxt(filename, dtype=None,comments="#",names=['lat','lon','pixsize','val'])
        ierror=0
    except IOError:
        print "File ("+filename+") does not exist in the current directory - "+currentdir
        ierror=0
    
    return ierror,epixarr
    
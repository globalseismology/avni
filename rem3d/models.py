#!/usr/bin/env python
"""
This script/module contains routines that are used to analyze Earth models and files that
contain them.
"""

#####################  IMPORT STANDARD MODULES   ######################################   

import os
import numpy as np #for numerical analysis

#####################

def readepixfile(filename):
    """Read .epix file format."""

    currentdir=os.getcwd()
    try: 
        f = open(filename, 'r')
    except IOError:
        print "File ("+filename+") does not exist in the current directory - "+currentdir
        sys.exit(2)
            
    epixarr=np.genfromtxt(filename, dtype=None,comments="#",names=['lat','lon','pixsize','val'])

    return epixarr
    
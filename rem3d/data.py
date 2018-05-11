#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format. 
Author: Raj Moulik, 2018"""

#####################  IMPORT STANDARD MODULES   ######################################   

import sys,os
import argparse #parsing arguments
import glob # pattern matching for finding files
import numpy as np #for numerical analysis
from scipy import stats #for stats analysis
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output and compare times
import pdb    #for the debugger pdb.set_trace()
from math import pi
import fortranformat as ff #reading/writing fortran formatted text
import requests
import platform
from datetime import datetime

####################### IMPORT REM3D LIBRARIES  #######################################

from . import constants
from . import tools

#######################################################################################

def creation_date(path_to_file):
    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    Modified from 
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    to use datetime
    """
    if platform.system() == 'Windows':
        return datetime.fromtimestamp(os.path.getctime(path_to_file))
    else:
        stat = os.stat(path_to_file)
        try:
            return datetime.fromtimestamp(stat.st_birthtime)
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return datetime.fromtimestamp(stat.st_mtime)
            
def update_file(installdir,file):
    """
    Does the url contain a downloadable resource that is newer
    """
    localfile = tools.get_filedir(checkwrite=True)+'/'+file
    url = constants.downloadpage + '/'+file
    h = requests.head(url, allow_redirects=True)
    if h.status_code == 404:
        print "Warning: Unknown status code ("+str(h.status_code)+") while quering "+file
        download = False
    elif h.status_code == 200:
        header = h.headers
        lmd = header.get('Last-Modified')  # Check when the file was modified
        server_data = datetime.strptime(lmd, '%a, %d %b %Y %H:%M:%S %Z')
        if os.path.isfile(localfile): # if a local file already exists
            local_data = creation_date(localfile)
            if (server_data-local_data).total_seconds() > 0: download = True # Download if server has newer file. 
        else:
            download = True
    else:
        print "Warning: Unknown status code ("+str(h.status_code)+") while quering "+file
        download = False
        
    if download:
        print ".... Downloading "+file+" from "+url
        r = requests.get(url, allow_redirects=True)
        open(localfile, 'wb').write(r.content)
    return 

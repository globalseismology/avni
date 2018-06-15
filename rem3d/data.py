#!/usr/bin/env python

"""
This script/module contains routines that are used to analyze data and files that
contain them.
"""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import os
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
    Modified from
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    to use datetime
    
    Parameters
    ----------

    path_to_file : full path to a file 
    
    Return
    ----------
    
    datetime stamp in UTC as REM3D server stores datetime in UTC

    """
    if platform.system() == 'Windows':
        return datetime.utcfromtimestamp(os.path.getctime(path_to_file))
    else:
        stat = os.stat(path_to_file)
        try:
            return datetime.utcfromtimestamp(stat.st_birthtime)
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            return datetime.utcfromtimestamp(stat.st_mtime)

def update_file(file):
    """
    If the REM3D url contain a downloadable resource that is newer, download it locally.
    """
    localfile = tools.get_filedir(checkwrite=True,makedir=True)+'/'+file
    url = constants.downloadpage + '/'+file
    h = requests.head(url, allow_redirects=True)
    download=False
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
        print ".... Downloading "+file+" from "+url+" to "+localfile
        r = requests.get(url, allow_redirects=True)
        open(localfile, 'wb').write(r.content)

    return

    
#!/usr/bin/env python

"""
This script/module contains routines that are used to analyze data and files that
contain them.
"""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

import os
import requests
import platform
from datetime import datetime
import calendar
import warnings

####################### IMPORT AVNI LIBRARIES  #######################################

from .. import constants
from .. import tools

#######################################################################################

def creation_date(path_to_file: str):
    """Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.

    Parameters
    ----------
    path_to_file : str
        full path to a file

    Returns
    -------
    datetime
        datetime stamp in UTC as AVNI server stores datetime in UTC

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2020.01.06 11.00
    """


    if platform.system() == 'Windows':
        return datetime.utcfromtimestamp(os.path.getctime(path_to_file))
    else:
        stat = os.stat(path_to_file)
        # We're probably on Linux. No easy way to get creation dates here,
        # so we'll settle for when its content was last modified.
    return datetime.utcfromtimestamp(stat.st_mtime)

def update_file(file,folder = None, baseurl = None, subdirectory = None):
    """If the AVNI server contain a downloadable resource that is newer, download it locally.

    Parameters
    ----------
    file : str
        full path and name of the file to sync with AVNI servers
    folder : str
        folder where local files are store, by default as output from tools.get_filedir()
    baseurl: str
        public URL from where the public downloads can take place, by default as specified as `downloadpage` in `constants.py`
    subdirectory: str
        subdirectory inside the baseurl where the file should be synced

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.01.06 11.00
    """

    # Get the correct default paths
    if folder is None: folder = tools.get_filedir()
    if baseurl is None: baseurl = constants.downloadpage
    if subdirectory is not None: baseurl = baseurl+'/'+subdirectory

    # Figure out internal-external path links
    localfile = os.path.join(folder,file)
    url = os.path.join(baseurl,file)

    # Perform http request
    h = requests.head(url, allow_redirects=True)
    download=False
    success = False
    if h.status_code == 404:
        warnings.warn("Warning: File not found with status code ("+str(h.status_code)+") while querying "+url)
    elif h.status_code == 200:
        header = h.headers

        # Check when the file was modified
        lmd = header.get('Last-Modified')
        server_data = datetime.strptime(lmd, '%a, %d %b %Y %H:%M:%S %Z')

        # Check if a local file already exists
        if os.path.isfile(localfile):
            local_data = creation_date(localfile)

            # Download if server has newer file
            if (server_data-local_data).total_seconds() > 0:
                download = True
            else:
                success = True
        else:
            download = True
    else:
        warnings.warn("Warning: Unknown status code ("+str(h.status_code)+") while querying "+file)

    # download if needed
    if download:
        print(".... Downloading "+file+" from AVNI server to "+localfile)
        r = requests.get(url, allow_redirects=True)
        open(localfile, 'wb').write(r.content)

        # calendar assumes tuple in UTC
        utime = calendar.timegm(server_data.timetuple())

        # set the time to server time but taking UTC into account above
        os.utime(localfile, (utime, utime))
        print(".... Download completed.")
        success = True
    return localfile,success


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
import pandas as pd
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
        print ".... Downloading "+file+" from "+url
        r = requests.get(url, allow_redirects=True)
        open(localfile, 'wb').write(r.content)

    return

def get_info_datafile(filename,extension='.interp.'):
    "Read the information from a file name while excluding file extension"
    typestr = filename.split('/')[-1].split(extension)[0]
    group, overtone, wavetype, period =typestr.split('.')
    return group, overtone, wavetype, period 
    
    
##################################   Surface waves #######################################

def read_SWhitcount(hitfile):
    """Read the the surface wave hit count file from sw_hitcount. The hit counts are normalized by area
     i.e. the values in the pixels with smallest ares are upweighted. """
    hit_array = np.genfromtxt(tools.get_fullpath(hitfile),dtype = None,names = ['lat','lon','val'],comments = "#")
    # grid spacing assuming a even grid
    grid_spacing = min(max(np.ediff1d(hit_array['lat'])),max(np.ediff1d(hit_array['lon'])))
    return hit_array,grid_spacing

def readREM3DSWformat(file,use_pandas=True):
    """Reads the REM3D format for analysis and plotting"""

    if (not os.path.isfile(file)): sys.exit("Filename ("+file+") does not exist")
    namelist=['overtone','peri','typeiorb','cmtname','eplat','eplon','cmtdep', \
              'stat','stlat','stlon','distkm','refphase','delobsphase', \
              'delerrphase','delpredphase']
    formatlist=['i','f8','a2','a15','f8','f8',  \
                'f8','a8','f8','f8','f8','f8',  \
                'f8','f8','f8']

    dtype = dict(names = namelist, formats=formatlist)

    # checks for CITE and SHORTCITE comments
    comments=[]
    for line in open(file):
        if line.startswith("#"): comments.append(line)
    reference='None'
    for line in open(file):
        if line.startswith("#CITE:"):
            reference=line[6:].rstrip('\n')
            break
    if reference=='None': print "Reference (e.g. #CITE: Ekstrom, 2011) should be defined as a comment in "+file ; sys.exit(2)
    shortref='None'
    for line in open(file):
        if line.startswith("#SHORTCITE:"):
            shortref=line[11:].rstrip('\n')
            break
    if shortref=='None': print "Short reference name (e.g. #SHORTCITE: GDM52) should be defined as a comment in "+file ; sys.exit(2)

    if use_pandas:
        # reads into a pandas dataframe (faster)
        SWdata=pd.read_csv(file, sep='\s+',names=namelist,comment='#')
    else:
        # reads into numpy structured array (slower)
        SWdata = np.genfromtxt(file, dtype=dtype,names=namelist,comments="#")

    return SWdata,comments,reference,shortref


def writeREM3DSWformat(filename,SWdata,comments):
    """Writes the REM3D format for analysis and plotting"""

    header_line = ff.FortranRecordWriter('(i2,1x,f6.1,1x,a2,1x,a15,1x,f8.3,1x,f9.3,1x, \
                                           f6.1,1x,a8,1x,f8.3,1x,f9.3,1x,f9.1,1x,f10.3,1x, \
                                           f10.3,1x,f10.3,1x,f10.3)')
    printstr = list(comments)
    for ii in np.arange(len(SWdata['overtone'])):
        arow=header_line.write([SWdata['overtone'][ii],SWdata['peri'][ii],SWdata['typeiorb'][ii],SWdata['cmtname'][ii].ljust(15),SWdata['eplat'][ii],SWdata['eplon'][ii],SWdata['cmtdep'][ii],SWdata['stat'][ii].ljust(8),SWdata['stlat'][ii],SWdata['stlon'][ii],SWdata['distkm'][ii],SWdata['refphase'][ii],SWdata['delobsphase'][ii],SWdata['delerrphase'][ii],SWdata['delpredphase'][ii]])
        printstr.append(arow+'\n')
    f=open(filename,'w')
    f.writelines(printstr)
    f.close()
    print "....written "+str(len(SWdata))+" observations to "+filename
    return

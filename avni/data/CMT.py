#!/usr/bin/env python

#######################################################################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

#####################  IMPORT STANDARD MODULES   ######################################

import os
import sys
import glob
import numpy as np
from datetime import datetime

if sys.version_info[0] >= 3: unicode = str

####################### IMPORT AVNI LIBRARIES  #######################################

from .common import update_file
from ..f2py import loadnbn2memory, getcmtbyname

####################       I/O ROUTINES     ######################################
def update_gcmt_nbn(standard='allorder.nbn',quick='qcmt.nbn'):
    """Updates the NBN file containing moment tensors from the Global CMT Project.

    Input parameters:
    ----------------

    standard :  file used in the standard catalogs of published CMTs

    quick :  file containing quick and unpublished CMTs from recent earthquakes

    """
    file_standard_gcmt = update_file(standard)
    file_quick_gcmt = update_file(quick)
    return file_standard_gcmt,file_quick_gcmt

def load_gcmt_nbn(choice=1, standard='allorder.nbn',quick='qcmt.nbn'):
    """Updates the NBN file containing moment tensors from the Global CMT Project.

    Input parameters:
    ----------------

    choice: 1 for both standard/quick (default), 2 for standard, 3 for quick

    standard :  file used in the standard catalogs of published CMTs

    quick :  file containing quick and unpublished CMTs from recent earthquakes

    """
    file_standard_gcmt,file_quick_gcmt = update_gcmt_nbn(standard=standard,quick=quick)
    nread,ierror = loadnbn2memory(choice,file_standard_gcmt,file_quick_gcmt)
    if ierror != 0: raise IOError('Could not read nbn files from the Global CMT catalog')
    return nread

def get_gcmt_info(cmtname, prefixes=['J','C']):
    """Updates the NBN file containing moment tensors from the Global CMT Project.

    Input parameters:
    ----------------

    cmtname :  CMTNAME of the earthquake in the Global CMT catalog.

    prefix : Some prefixes that are also checked if cmtname is not found.

    Output:
    ----------------

    time : time of earthquake in datetime format, precision to microsecond

    elat, elon, edep : latitude, longitude and depth of earthquake

    Mw : Moment magnitude of earthquake

    ievt : event index

    """
    ievt,iyear,month,iday,ihour,minute,fsec, \
        elat,elon,edep,Mw,ierror = getcmtbyname(cmtname)

    icount=0
    while ierror != 0 and icount < len(prefixes):
        prefix = prefixes[icount]
        ievt,iyear,month,iday,ihour,minute,fsec, \
            elat,elon,edep,Mw,ierror = getcmtbyname(prefix+cmtname)
        icount += 1

    if ierror != 0:
        raise IOError('No earthquakes found in Global CMT catalog in getcmtdate for '+cmtname)
    else:
        time = datetime(iyear,month,iday,ihour,minute, \
            int(fsec),microsecond=int((fsec-int(fsec))*1000000))
    return time,elat,elon,edep,Mw,ievt
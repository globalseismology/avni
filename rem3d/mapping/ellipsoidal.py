#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function

import multiprocessing
from joblib import Parallel, delayed
import pdb
############################### PLOTTING ROUTINES ################################
from ..tools.trigd import atand,tand
from ..f2py import ddelazgc # geolib library from NSW
from .. import constants
###############################

def get_distaz(eplat,eplon,stlat,stlon,num_cores=1):
    """Get the distance and azimuths from positions in geographic coordinates"""

    geoco = constants.geoco.magnitude
    if isinstance(eplat,list): # if the input is a list loop
        delta=[];azep=[];azst=[]
        # Standard checks on number of cores
        avail_cores = multiprocessing.cpu_count()
        if num_cores > avail_cores:
            raise ValueError("Number of cores requested ("+str(num_cores)+") is higher than available ("+str(avail_cores)+")")
        # Crate list of tuples of job arguments and pass to parallel using a helper routine
        job_args = [(atand(geoco*tand(item_lat)),eplon[jj],atand(geoco*tand(stlat[jj])),stlon[jj]) for jj, item_lat in enumerate(eplat)]
        temp=Parallel(n_jobs=num_cores)(delayed(delazgc_helper)(ii) for ii in job_args)
        for il in temp: delta.append(il[0]);azep.append(il[1]);azst.append(il[2])

    elif isinstance(eplat,float):
        elat=atand(geoco*tand(eplat))
        elon=eplon
        slat=atand(geoco*tand(stlat))
        slon=stlon
        delta,azep,azst = ddelazgc(elat,elon,slat,slon)
    else:
        raise ValueError("get_distaz only takes list or floats")

    return delta,azep,azst

def delazgc_helper(args):
    return ddelazgc(*args)
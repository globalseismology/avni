#!/usr/bin/env python

"""
This script/module contains routines that are used to analyze Earth models and files that
contain them.
"""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import os,sys
import numpy as np #for numerical analysis

####################### IMPORT REM3D LIBRARIES  #######################################

from . import tools   

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
    except IOError:
        raise IOError("File (",filename,") does not exist in the current directory - ",currentdir)
    
    return epixarr

def read3dmodelfile(modelfile,maxkern=300,maxcoeff=6000):
    """
    Reads a standard 3D model file. maxkern is the maximum number of radial kernels
    and maxcoeff is the maximum number of corresponding lateral basis functions.
    resolution and realization are the indices for the resolution level
    and the realization from a model ensemble (usually 0 if a single file)
    """    
    if (not os.path.isfile(modelfile)): raise IOError("Filename ("+modelfile+") does not exist")

    desckern=np.zeros(maxkern, dtype='a40')
    ihorpar=np.zeros(maxkern, dtype=int)
    coef=np.zeros([maxcoeff,maxkern], dtype=float)
    with open(modelfile) as f: lines = f.readlines()
    ii=0
    while ii < len(lines):
        line=lines[ii]; ii=ii+1
        if line.startswith("REFERENCE MODEL:"): refmodel=line[16:].rstrip('\n').strip(' ')
        if line.startswith("KERNEL SET:"): kerstr=line[12:].rstrip('\n').strip(' ')
        if line.startswith("RADIAL STRUCTURE KERNELS:"): nmodkern = int(line[26:].rstrip('\n'))
        if line.startswith("DESC") and line[8] == ':': 
            idummy=int(line[4:8])
            substr=line[9:len(line.rstrip('\n'))]
            ifst,ilst=tools.firstnonspaceindex(substr)
            desckern[idummy-1]=substr[ifst:ilst]
        if line.startswith("HORIZONTAL PARAMETERIZATIONS:"): 
            nhorpar = int(line[29:].rstrip('\n'))
            hsplfile=np.zeros(nhorpar, dtype='a40')
            typehpar=np.zeros(nhorpar, dtype='a40')
            ityphpar=np.zeros(nhorpar, dtype=int)
            lmaxhor=np.zeros(nhorpar, dtype=int)
            ncoefhor=np.zeros(nhorpar, dtype=int)
            ixlspl=np.zeros([maxcoeff,nhorpar], dtype=int)
            xlaspl=np.zeros([maxcoeff,nhorpar], dtype=float)
            xlospl=np.zeros([maxcoeff,nhorpar], dtype=float)
            xraspl=np.zeros([maxcoeff,nhorpar], dtype=float)
        if line.startswith("HPAR") and line[8] == ':':     
            idummy=int(line[4:8])
            substr=line[9:len(line.rstrip('\n'))]
            ifst,ilst=tools.firstnonspaceindex(substr)
            if substr[ifst:ifst+20] == 'SPHERICAL HARMONICS,':
                lmax = int(substr[21:].rstrip('\n'))
                ityphpar[idummy-1]=1
                typehpar[idummy-1]='SPHERICAL HARMONICS'
                lmaxhor[idummy-1]=lmax
                ncoefhor[idummy-1]=(lmax+1)**2
            elif substr[ifst:ifst+18] == 'SPHERICAL SPLINES,':
                ifst1=ifst+18
                ifst=len(substr)
                ilst=len(substr)
                while substr[ifst-1:ifst] != ',': ifst=ifst-1
                ncoef=int(substr[ifst+1:ilst].rstrip('\n'))
                substr=substr[ifst1:ifst-1]
                ifst1,ilst=tools.firstnonspaceindex(substr)
                hsplfile[idummy-1]=substr[ifst1:ilst]
                ityphpar[idummy-1]=2
                typehpar[idummy-1]='SPHERICAL SPLINES'
                lmaxhor[idummy-1]=0
                ncoefhor[idummy-1]=ncoef
                for jj in range(ncoef):
                    arr=lines[ii].rstrip('\n').split(); ii=ii+1
                    ixlspl[jj,idummy-1]=int(arr[0]); xlaspl[jj,idummy-1]=arr[1]
                    xlospl[jj,idummy-1]=arr[2]; xraspl[jj,idummy-1]=arr[3]
        if line.startswith("STRU") and line[8] == ':':
            idummy=int(line[4:8])
            ihor=int(line[9:].rstrip('\n'))
            ihorpar[idummy-1]=ihor
            ncoef=ncoefhor[ihor-1]
            for jj in range(ncoef/6):
                arr=lines[ii].rstrip('\n').split(); ii=ii+1
                coef[jj*6:(jj+1)*6,idummy-1]=[float(i) for i in arr]
            remain = ncoef % 6    
            arr=lines[ii].rstrip('\n').split(); ii=ii+1
            coef[(jj+1)*6:(jj+1)*6+remain,idummy-1]=[float(i) for i in arr]

    # Store the variables
    numvar=0; varstr=np.zeros(nmodkern, dtype='a40')
    ivarkern=np.zeros(nmodkern)
    for ii in np.arange(nmodkern):
        string=desckern[ii]
        #pdb.set_trace()
        for kk in np.arange(numvar):
            if varstr[kk] == string[:string.index(',')]:
                ivarkern[ii]=kk+1
        if ivarkern[ii] == 0:
            numvar=numvar+1
            varstr[numvar-1] = string[:string.index(',')]
            ivarkern[ii]=numvar

    # Save the relevant portions
    desckern = desckern[:nmodkern]
    ihorpar = ihorpar[:nmodkern]
    varstr = varstr[:numvar]
    coef = coef[:max(ncoefhor),:nmodkern].transpose() # to get it in a kernel * coeff format
    ixlspl = ixlspl[:max(ncoefhor),:].transpose() 
    xlaspl = xlaspl[:max(ncoefhor),:].transpose() 
    xlospl = xlospl[:max(ncoefhor),:].transpose() 
    xraspl = xraspl[:max(ncoefhor),:].transpose() 
    
    # Store in a dictionary
    model3d = {}
    model3d['refmodel']=refmodel; model3d['kerstr']=kerstr;model3d['nmodkern']=nmodkern
    model3d['desckern']=desckern; model3d['nhorpar']=nhorpar;model3d['hsplfile']=hsplfile
    model3d['ityphpar']=ityphpar; model3d['typehpar']=typehpar; model3d['lmaxhor']=lmaxhor
    model3d['ncoefhor']=ncoefhor; model3d['ixlspl']=ixlspl; model3d['xlaspl']=xlaspl
    model3d['xlospl']=xlospl; model3d['xraspl']=xraspl; model3d['ihorpar']=ihorpar
    model3d['ivarkern']=ivarkern; model3d['numvar']=numvar
    model3d['varstr']=varstr; model3d['coef']=coef
    return model3d
    

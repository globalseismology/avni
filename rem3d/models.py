#!/usr/bin/env python

"""
This script/module contains routines that are used to analyze Earth models and files that
contain them.
"""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import numpy as np #for numerical analysis
import pdb    #for the debugger pdb.set_trace()
from math import pi

import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases
import matplotlib.pyplot as plt
from scipy import sparse
import re
from copy import copy, deepcopy
import struct
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

#####################
# 3D model class
class model3d(object):
    '''
    A class for 3D reference Earth models used in tomography
    '''

    def __init__(self,file=None,num_resolution=1,num_realization=1):
        self.metadata ={}
        self.data = {}
        for ii in np.arange(num_resolution): 
            self.metadata['resolution_'+str(ii)] = None
            self.data['resolution_'+str(ii)] = {}
            for jj in np.arange(num_realization):
                self.data['resolution_'+str(ii)]['realization_'+str(jj)] = {'name':None,'coef':None}
        self.name = None
        self.type = None
        self.refmodel = None
        self.description = None
        if file is not None: self.readfile(file)
    
    def __str__(self):
        if self.name is not None:
            output = "%s is a three-dimensional model ensemble with %s resolution levels, %s realizations of %s type and %s as the reference model" % (self.name,len(self.metadata), len(self.data['resolution_0']),self.type,self.refmodel)
        else:
            output = "No three-dimensional model has been read into this model3d instance yet"
        return output
        
                
    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

    def readfile(self,modelfile,resolution=0,realization=0,maxkern=300,maxcoeff=6000):
        """
        Reads a standard 3D model file. maxkern is the maximum number of radial kernels
        and maxcoeff is the maximum number of corresponding lateral basis functions.
        resolution and realization are the indices for the resolution level
        and the realization from a model ensemble (usually 0 if a single file)
        """    
        if (not os.path.isfile(modelfile)): raise IOError("Filename ("+modelfile+") does not exist")
    
        # read mean model   
        model=read3dmodelfile(modelfile)
    
        # Write to a dictionary
        metadata = {}
        metadata['kerstr']=model['kerstr']; metadata['nmodkern']=model['nmodkern']
        metadata['desckern']=model['desckern']; metadata['nhorpar']=model['nhorpar']; metadata['hsplfile']=model['hsplfile']
        metadata['typehpar']=model['typehpar']; metadata['ityphpar']=model['ityphpar']; metadata['lmaxhor']=model['lmaxhor']
        metadata['ixlspl']=model['ixlspl']; metadata['xlaspl']=model['xlaspl']; metadata['xlospl']=model['xlospl']
        metadata['xraspl']=model['xraspl']; metadata['ihorpar']=model['ihorpar']; metadata['ivarkern']=model['ivarkern']
        metadata['numvar']=model['numvar']; metadata['varstr']=model['varstr']
        metadata['ncoefhor']=model['ncoefhor']
        
        self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['name'] = ntpath.basename(modelfile)
        self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['coef'] = sparse.csr_matrix(model['coef'])

        self.metadata['resolution_'+str(resolution)] = metadata
        self.refmodel = model['refmodel']
        self.name = ntpath.basename(modelfile)
        self.description = "Read from "+modelfile
        self.type = 'rem3d'
        
        return 

    def coeff2modelarr(self,resolution=[0],realization=0):
        """
        Convert the coeffiecient matrix from the file to a sparse model array. Add
        up all the resolutions if a list is provided.
        
        realization : index of the set of coefficients in an ensemble. Default is 0
        as there is only one set of coefficients when read from a model file.
        
        resolution : list of resolutions to include the the modelarray
        """
        if self.type != 'rem3d': raise NotImplementedError('model format ',self.type,' is not currently implemented in reference1D.coeff2modelarr')
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        
        if isinstance(resolution, (list,tuple,np.ndarray)):
            resolution = np.asarray(resolution)
        else:
            raise TypeError('resolution must be function or array, not %s' % type(resolution))
            
        # Loop over resolution levels, adding coefficients 
        for ir in range(len(resolution)): 
            try:
                coefficients =self.data['resolution_'+str(resolution[ir])]['realization_'+str(realization)]['coef'].toarray()
            # Loop over all kernel basis
                for ii in range(len(self.metadata['resolution_'+str(resolution[ir])]['ihorpar'])): 
                    # number of coefficients for this radial kernel
                    ncoef = self.metadata['resolution_'+str(resolution[ir])]['ncoefhor'][self.metadata['resolution_'+str(resolution[ir])]['ihorpar'][ii]-1]
                    # first radial kernel and first tesselation level
                    if ii == 0 and ir == 0: 
                        modelarr=coefficients[ii,:ncoef]
                    else:
                        modelarr=np.append(modelarr,coefficients[ii,:ncoef]) 
            except AttributeError: # 
                raise ValueError('resolution '+str(resolution[ir])+' and realization '+str(realization)+' not filled up yet.')
        modelarr = sparse.csr_matrix(modelarr) # Convert to sparse matrix
        modelarr = modelarr.T # take the transpose to make dot products easy 
        return modelarr

    

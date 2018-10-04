#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format. 
Author: Raj Moulik, 2018"""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import argparse #parsing arguments
import glob # pattern matching for finding files
import numpy as np #for numerical analysis
from datetime import date  #to give a timestamp to output and compare times
import pdb    #for the debugger pdb.set_trace()
from math import pi
import fortranformat as ff #reading/writing fortran formatted text
from future.utils import native_str
from six import string_types # to check if variable is string using isinstance
from numpy.lib.recfunctions import append_fields # for appending named columns to named numpy arrays
from scipy.interpolate import griddata
import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases
import matplotlib.pyplot as plt
from matplotlib import gridspec # Relative size of subplots
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from scipy import sparse
from configobj import ConfigObj
import re
from copy import copy, deepcopy
import struct
import h5py

####################### IMPORT REM3D LIBRARIES  #######################################
from . import tools   
from . import plots 
from . import constants
from . import io                
#######################################################################################
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

    desckern=np.zeros(maxkern, dtype='U40')
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
            hsplfile=np.zeros(nhorpar, dtype='U40')
            typehpar=np.zeros(nhorpar, dtype='U40')
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
            for jj in range(int(ncoef/6)):
                arr=lines[ii].rstrip('\n').split(); ii=ii+1
                coef[jj*6:(jj+1)*6,idummy-1]=[float(i) for i in arr]
            remain = ncoef % 6    
            arr=lines[ii].rstrip('\n').split(); ii=ii+1
            coef[(jj+1)*6:(jj+1)*6+remain,idummy-1]=[float(i) for i in arr]

    # Store the variables
    numvar=0; varstr=np.zeros(nmodkern, dtype='U40')
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
    
def readensemble(filename):
    """Read .npz file containing a collection of models.

    Parameters
    ----------

    filename : An ensemble file

    """


    return 

def readprojections(type='radial'):
    """Read .npz file containing a collection of models.

    Parameters
    ----------

    filename : An file containing all projection matrices

    """


    return 
#####################
# Vertical basis parameter class that defines an unique combination of functions, their radial parameterization and any scaling
# 3D model class
class radial_basis(object):
    '''
    A class for radial bases that defines a unique combination of parameters,
    their radial parameterization and any scaling that is used.
    '''
    def __init__(self):
        self.proj = {'data': None,'to_name':None,'to_type':'boxcar','to_attributes':{}}
        self.metadata = {'from_name':None,'from_type':None,'from_attributes':{}}
    
    def readprojfile(self,projverfile):
        """
        Read projection matrix for going between vertical bases from a file
        """    

    def eval_radial(self,depth):
        """
        Evaluate radial basis at a depth interpolated from existing projection matrices.
        """    
    
    def project_boxdepth(self,depth_range):
        """
        Project from current vertical basis to a vertical boxcar basis
        depth_range is named numpy array of top and bottom depths
        """    
        
        
#####################
# Horizontal basis parameter class that defines an unique combination of parameters, their radial parameterization and any scaling
# 3D model class
class lateral_basis(object):
    '''
    A class for radial bases that defines a unique combination of parameters,
    their radial parameterization and any scaling that is used.
    '''
    def __init__(self, name, type):
        """
        types : 'epix','ylm','sh','wavelet','slepians'
        """
        self.proj = {}
        self.metadata = {'name':name,'type':type,'attributes':{}}

    def addtypes(self, names, types):
        """
        types = ['ylm','sh','wavelet','slepians']
        """
        # check the type 
        if isinstance(types,string_types): types = np.array(types) 
        if len(names) !=len(types): raise ValueError("len(names) !=len(types)")

        for ii in np.arange(len(types)):
            # if does not exist already`or not the same as self type
            if types[ii] not in self.proj.keys() and types[ii] != self.metadata['type']:         
                self.proj[types[ii]] = {}
                self.proj[types[ii]][to_name[ii]] = {'data':None,'attributes':{}}        
    
    def readprojfile(self,infile):
        """
        Reads a projection matrix file that evaluates the radial bases at various depths.
        """
        
        if (not os.path.isfile(infile)): raise IOError("Filename ("+infile+") does not exist. Use shell script print_projmatrix to create it.")
        nbytes = os.path.getsize(infile)

        cc = 0 #initialize byte counter
        ifswp = '' # Assuem that byte order is not be swapped unless elat is absurdly high

        with open(infile, "rb") as f:
            # preliminary metadata
            indata = f.read(4); cc = cc+4 # try to read iflag
            iflag = struct.unpack(ifswp+'i',indata)[0] # Read flag
            if iflag != 1: 
                ifswp = '!' # swap endianness from now on
                iflag = struct.unpack(ifswp+'i',indata)[0]
                if iflag != 1: sys.exit("Error: iflag != 1")
            self.metadata['from_type'] = struct.unpack('20s',f.read(20))[0].strip(); cc = cc+20
            pdb.set_trace()
            self.proj['to_type'] = struct.unpack('20s',f.read(20))[0].strip(); cc = cc+20
            #ndp = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
            #npx = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
            #nhorcum = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4

        
    def eval_lateral(self,lat,lon):
        """
        Evaluate radial basis at a depth interpolated from existing projection matrices.
        """    

    def project_lateral(self,lat,lon):
        """
        Project from current horizontal basis to another orthogonal basis 
        and return the coefficients.
        """    

#####################
# Kernel set class that stores the list 
# 3D model class
class radial_basis(object):
    '''
    A class for radial bases that defines a unique combination of parameters,
    their radial parameterization and any scaling that is used.
    '''
    def __init__(self):
        self.data = None
        self.metadata = {}
        self.metadata['name'] = None
        self.metadata['type'] = None
        self.metadata['depths'] = None
    
    def readprojfile(self,projfile):
        """
        Reads a projection matrix file that evaluates the radial bases at various depths.
        """    

    def eval_radial(self,depth):
        """
        Reads a projection matrix file that evaluates the radial bases at various depths.
        """    

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
      
    def add_realization(self,resolution=None,name=None,coef=None):
        """
        Added a set of realizations to the object. resolution is the tesselation level at
        which the instance is update with name and coeff (default: None).
        """      
        if resolution==None:
            resolution = len(self.data)
            for ii in np.arange(resolution):
                realization = len(self.data['resolution_'+str(ii)]) 
                self.data['resolution_'+str(ii)]['realization_'+str(realization)] = {'name':name,'coef':coef}
        else:
            if isinstance(resolution, int) :
                realization = len(self.data['resolution_'+str(resolution)])
                self.data['resolution_'+str(resolution)]['realization_'+str(realization)] = {'name':name,'coef':coef}    
            else:
                raise ValueError('Invalid resolution in add_realization')

    def add_resolution(self):
        """
        Added a resolution level to the object. num_realization is the number 
        of realizations of model coefficients within that object.
        """      
        num_resolution = len(self.data)
        self.metadata['resolution_'+str(num_resolution)] = None
        self.data['resolution_'+str(num_resolution)] = {}
        self.data['resolution_'+str(num_resolution)]['realization_0'] = {'name':None,'coef':None}

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
        
    def getradialattributes(self,parserfile='attributes.ini'):
        """
        Reads configuration file and get the basis attributes like knot depth for each 
        parameter in modelarray list of coefficients. parser is the output from 
        parser = ConfigObj(config) where config is the configuration *.ini file.
        """    
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        
        # Read configuration file to a configuration parser. This is to make this available on the fly
        filepath = tools.get_configdir(checkwrite=True) + '/'+parserfile
        if (not os.path.isfile(filepath)): raise IOError("Configuration file ("+filepath+") does not exist")
        parser = ConfigObj(filepath)
    
        for ll in np.arange(len(self.metadata)): # Loop over every resolution

            # Read the kernel set from the model3d dictionary 
            kerstr=self.metadata['resolution_'+str(ll)]['kerstr']

            # Read the basis
            temp = parser['Kernel_Set'][kerstr]['radial_knots']; radial_knots = []
            for ii in range(len(temp)): 
                temp2 = [float(i) for i in temp[ii].split(',')]
                radial_knots.append(temp2)

            # Clustering configuration
            radial_type = parser['Kernel_Set'][kerstr]['radial_type']            
        
            # Loop over all kernel basis
            knot_depth=[]
            for ii in range(len(self.metadata['resolution_'+str(ll)]['ihorpar'])): 
                ifound=0
                for jj in range(len(radial_type)):
                    if re.search(radial_type[jj],self.metadata['resolution_'+str(ll)]['desckern'][ii]):
                        index = int(self.metadata['resolution_'+str(ll)]['desckern'][ii].split(',')[-1]) - 1
                        knot_depth.append(radial_knots[jj][index]); ifound=1
                if ifound != 1: 
                    print("Warning: Did not find radial kernel knot depth in getradialattributes for "+self.metadata['resolution_'+str(ll)]['desckern'][ii]+". Setting to NaN to denote ignorance in clustering")
                    knot_depth.append(np.nan)
            # Stor in relevant variables
            self.metadata['resolution_'+str(ll)]['knot_depth']=np.array(knot_depth)
            self.metadata['resolution_'+str(ll)]['description']=parser['Kernel_Set'][kerstr]['description']
            self.metadata['resolution_'+str(ll)]['scaling']=parser['Kernel_Set'][kerstr]['scaling']

        return 
            
    def readprojmatrix(self,lateral_basis):
        """
        Reads Projection matrix created by plot_3dmod_pm. 
        lateral_basis can be M362 or pixel1
        """    
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")

        #read all the bytes to indata
        outfile = self.name+'.'+lateral_basis+'.proj.bin.h5'
        infile = self.name+'.'+lateral_basis+'.proj.bin'
    
        xlat=[];xlon=[];area=[];deptharr=[]
        cc = 0 #initialize byte counter
        ifswp = '' # Assuem that byte order is not be swapped unless elat is absurdly high

        if (not os.path.isfile(outfile)): 
            if (not os.path.isfile(infile)): raise IOError("Filename ("+infile+") does not exist. Use shell script plot_3dmod_pm64 to create it.")
            nbytes = os.path.getsize(infile)
            print("....writing "+outfile)
            with open(infile, "rb") as f:
                # preliminary metadata
                indata = f.read(4); cc = cc+4 # try to read iflag
                iflag = struct.unpack(ifswp+'i',indata)[0] # Read flag
                if iflag != 1: 
                    ifswp = '!' # swap endianness from now on
                    iflag = struct.unpack(ifswp+'i',indata)[0]
                    if iflag != 1: sys.exit("Error: iflag != 1")
                ndp = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
                npx = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
                nhorcum = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
                for jj in np.arange(npx):
                    xlat.append(struct.unpack(ifswp+'f',f.read(4))[0]); cc = cc+4
                    xlontemp = struct.unpack(ifswp+'f',f.read(4))[0]; cc = cc+4
                    # Change from -180 to 180 limits to be compatible with model3d format
                    if xlontemp > 180: xlontemp = xlontemp - 360.
                    xlon.append(xlontemp)
                    area.append(struct.unpack(ifswp+'f',f.read(4))[0]); cc = cc+4
                # loop by deptharr
                projarr={};refstrarr=[];refvalarr={}
                for ii in np.arange(ndp): 
                    deptharr.append(struct.unpack(ifswp+'f',f.read(4))[0]); cc = cc+4
                    for jj in np.arange(npx):
                        if jj % (npx/5) == 0:
                            print("deptharr # "+str(ii+1)+" out of "+str(ndp)+", pixel # "+str(jj+1)+" out of "+str(npx))
                        neval = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
                        for kk in np.arange(neval):
                            refstr = struct.unpack('80s',f.read(80))[0].strip(); cc = cc+80
                            refval=struct.unpack(ifswp+'f',f.read(4))[0]; cc = cc+4
                            if ii==0 and jj==0: refstrarr.append(refstr)
                            if jj==0: refvalarr[ii,kk]=refval
                            iloop=struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
                            coeirow=[];proja=[]
                            for ll in np.arange(iloop):
                                coeirow.append(struct.unpack(ifswp+'i',f.read(4))[0]); cc = cc+4
                                proja.append(struct.unpack(ifswp+'f',f.read(4))[0]); cc = cc+4
                            coeirow=np.array(coeirow)
                            proja=np.array(proja)
                            try:
                                projarr[ii,kk]=projarr[ii,kk]+sparse.csr_matrix( (proja,(np.ones_like(coeirow)*jj,coeirow-1)), shape=(npx,nhorcum))
                            except KeyError:
                                projarr[ii,kk]=sparse.csr_matrix( (proja,(np.ones_like(coeirow)*jj,coeirow-1)), shape=(npx,nhorcum))
            if cc != nbytes: sys.exit("Error: number of bytes read "+str(cc)+" do not match expected ones "+str(nbytes))
            deptharr=np.array(deptharr); refstrarr=np.array(refstrarr)
            # save grid
            namelist = ['xlat', 'xlon', 'area']
            namelist = [str(name) for name in namelist]
            formatlist=['f8','f8','f8']
            xlat=np.array(xlat); xlon=np.array(xlon); area=np.array(area)
            results = [xlat, xlon, area]        
            results = map(list, zip(*results))
            dtype = dict(names = namelist, formats=formatlist)     
            grid = np.array([tuple(x) for x in results],dtype=dtype)
            model = self.name
                        
            h5f = h5py.File(outfile,'w')
            h5f.attrs["infile"] = np.string_(ntpath.basename(infile))
            h5f.attrs["ndp"] = ndp
            h5f.attrs["npx"] = npx
            h5f.attrs["nhorcum"] = nhorcum
            h5f.attrs["neval"] = neval
            h5f.attrs["deptharr"] = deptharr
            h5f.attrs["refstrarr"] = refstrarr
            io.store_numpy_hdf(h5f,'grid',grid)
            h5f.attrs["model"] = np.string_(model)
            h5f.attrs["param"] = np.string_(lateral_basis)
            
            for ii in np.arange(len(deptharr)):
                for jj in np.arange(len(refstrarr)):
                    proj = h5f.create_group("projection/depth_"+str(ii)+ "/refstr_"+str(jj))
                    proj.attrs['refvalue']=refvalarr[ii,jj]
                    io.store_sparse_hdf(h5f,"projection/depth_"+str(ii)+ "/refstr_"+str(jj),projarr[ii,jj])
            h5f.close()   
        else:
            print("....reading "+outfile)
            h5f = h5py.File(outfile,'r')
            ndp=h5f.attrs['ndp']; npx=h5f.attrs['npx']; nhorcum=h5f.attrs['nhorcum']
            neval=h5f.attrs['neval']; deptharr=h5f.attrs['deptharr']
            refstrarr=h5f.attrs['refstrarr']; grid = io.load_numpy_hdf(h5f,'grid') 
            xlat=grid['xlat']; xlon=grid['xlon']; area=grid['area']
            model=h5f.attrs['model']; param=h5f.attrs['param']
            
            # Always read the following from hdf5 so projarr is list ofsparsehdf5-fast I/O
            projarr={};refvalarr={}
            for ii in np.arange(len(deptharr)):
                for jj in np.arange(len(refstrarr)):
                    proj = h5f["projection/depth_"+str(ii)+ "/refstr_"+str(jj)]
                    refvalarr[ii,jj] = proj.attrs['refvalue']
                    projarr[ii,jj] = io.load_sparse_hdf(h5f,"projection/depth_"+str(ii)+ "/refstr_"+str(jj))
            h5f.close()  
     
        # Write to a dictionary
        projection = {}
        projection['ndp']=ndp; projection['npx']=npx; projection['nhorcum']=nhorcum; 
        projection['neval']=neval; projection['deptharr']=deptharr; projection['refstrarr']=refstrarr
        projection['xlat']=xlat; projection['xlon']=xlon; projection['area']=area
        projection['refvalarr']=refvalarr; projection['projarr']=projarr       
        projection['model']=model; projection['param']=lateral_basis         
        return projection

    def getprojmatrix(self,to_name='epix',depths = [30.,50.]):
        """
        Get the projection matrix from a lateral basis to another and for particular depths  
        """    
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        
        # Get the radial projection file
        projfile,exists = tools.get_projections(type='radial')
        
        pdb.set_trace()
        
     
        # Write to a dictionary
        projection = {}
        projection['ndp']=ndp; projection['npx']=npx; projection['nhorcum']=nhorcum; 
        projection['neval']=neval; projection['deptharr']=deptharr; projection['refstrarr']=refstrarr
        projection['xlat']=xlat; projection['xlon']=xlon; projection['area']=area
        projection['refvalarr']=refvalarr; projection['projarr']=projarr       
        projection['model']=model; projection['param']=lateral_basis         
        return projection


    def getprojtimesmodel(self,projection,variable,depth,realization=0):
        """
        Projection matrix multiplied by model ensemble. Choses the nearest depth available for the projection.
        """    
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")

        modelarr = self.coeff2modelarr(realization=realization)
        # Select appropriate arrays from projection matrix, read from file
        projarr = projection['projarr']; refstrarr = projection['refstrarr']; deptharr = projection['deptharr']  
    
        varindex = np.where(refstrarr==variable)[0]
        absdiff=abs(deptharr-depth)
        # Find the nearest depth
        depindex = absdiff.argmin()
        if min(absdiff) > 0. :
            print ("No unique depth found in the projection matrix. Choosing the nearest available depth "+str(deptharr[depindex]))
        modelselect=projarr[depindex,varindex[0]]*modelarr    
        return modelselect,deptharr[depindex]

    def printsplinefiles(self):
        """
        Prints out the splines knots into a file.
    
        Parameters
        -----------
    
        model3d : The model object from read3dmodelfile
    
        """
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        for ii in np.arange(len(self.metadata)):
        
            ifindspline=np.where(self.metadata['resolution_'+str(ii)]['typehpar']=='SPHERICAL SPLINES')[0]
            for ifind in ifindspline:
                filename=self.metadata['resolution_'+str(ii)]['hsplfile'][ifind]
                arr = np.vstack([self.metadata['resolution_'+str(ii)]['xlospl'][ifind],self.metadata['resolution_'+str(ii)]['xlaspl'][ifind],self.metadata['resolution_'+str(ii)]['xraspl'][ifind]]).transpose()
                np.savetxt(filename,arr, fmt='%7.3f %7.3f %7.3f')   
                print(".... written "+filename)
        return


    def plotslices(self,lateral_basis='pixel1',dbs_path='~/dbs',x=0,percent_or_km='%',colormin = -6.,colormax=6.,depth=None):
        """
        Plots interactively a model slice of a variable at a given depth till an 
        invalid depth is input by the user
    
        Parameters
        ----------
        model3d : the model dictionary read by read3dmodelfile
        
        param : lateral parameterization dictionary read by readprojmatrix
        
        x,percent_or_km, colormin,colormax,depth : plotting options for jupyter 
                                                   instead of interactive input
        """    
        # Select appropriate modelarray m
        modelarr = self.coeff2modelarr() # convert to modelarr of d = G*m
    
        # Select appropriate arrays from projection matrix, read from file
        projection = self.readprojmatrix(lateral_basis)
        projarr = projection['projarr']; refstrarr = projection['refstrarr']; deptharr = projection['deptharr']  
        xlat = projection['xlat']; xlon = projection['xlon']
    
        # select models based on parameter and depth desired
        new_figure='y'  # flag for done
        while (new_figure =='y' or new_figure == 'Y'):
            plt.ion()
            fig=plt.figure() 
            try: 
                subplotstr = input("Provide rows and colums of subplots - default  is 1 1:")
                subploty,subplotx = int(subplotstr.split()[0]),int(subplotstr.split()[1])
            except (ValueError,IndexError,SyntaxError,EOFError):
                subploty = 1; subplotx=1

            flag=0  # flag for depth
            while (flag < subploty*subplotx):
                flag=flag+1
                try:
                    for ii in np.arange(len(refstrarr)): print(ii,refstrarr[ii])
                    try:
                        x = int(input("Select variable to plot - default is 0:"))
                    except (ValueError,EOFError):
                        x = x
                    try:
                        depth = float(input("Select depth - select any value for topography ["+str(min(deptharr))+"-"+str(max(deptharr))+"] :"))
                    except (ValueError,EOFError):
                        if depth is None:
                            depth = min(deptharr)
                        else:
                            depth = depth
                    if depth < min(deptharr) or depth > max(deptharr):
                        flag=flag-1
                    else:
                        modeleselect,_ = self.getprojtimesmodel(projection,variable=refstrarr[x],depth=depth)
                        try:
                            percent_or_km = input("Is the value in percentage or km or percentage -default is % [km/%]:")
                        except (EOFError):
                            percent_or_km = percent_or_km
                        if percent_or_km == '': percent_or_km = '%'                            
                        if percent_or_km == '%': 
                            modelarray = modeleselect.toarray().flatten()*100. # In Percent assuming elastic parameters
                        else:
                            modelarray = modeleselect.toarray().flatten()

                        # Get limits for colorbar
                        try: 
                            colorstr = input("Input two values for minimum and maximum values of colorbar - default is -6 6:")
                            colormin,colormax = float(colorstr.split()[0]),float(colorstr.split()[1])
                        except (ValueError,IndexError,EOFError):
                            colormin = colormin; colormax=colormax
            
                        # Plot the model
                        test = np.vstack((xlat,xlon,modelarray)).transpose()
                        dt = {'names':['lat', 'lon', 'val'], 'formats':[np.float, np.float, np.float]}
                        plotmodel = np.zeros(len(test), dtype=dt)
                        plotmodel['lat'] = test[:,0]; plotmodel['lon'] = test[:,1]; plotmodel['val'] = test[:,2]
                        ax=fig.add_subplot(subploty,subplotx,flag)
                        plots.globalmap(ax,plotmodel,colormin,colormax,dbs_path=dbs_path, colorlabel='Anomaly '+'(in '+percent_or_km+')' if percent_or_km!='' else 'Anomaly', grid=[30.,90.],gridwidth=0,projection='robin',lat_0=0, lon_0=150., colorpalette='rem3d',colorcontour=21)
                        ax.set_title(refstrarr[x]+' at '+str(depth)+' km.' if 'Topo' not in refstrarr[x] and 'topo' not in refstrarr[x] else refstrarr[x])
                        fig.canvas.draw() 
                except SyntaxError:
                    flag=flag-1
            try:
                new_figure = input("Another figure? y/n:")
            except (EOFError):
                new_figure = 'n'
    
        return 

#####################
# 1D model class

class reference1D(object):
    '''
    A class for 1D reference Earth models used in tomography
    '''

    def __init__(self):
        self.__nlayers__ = None
        self.data = None
        self.metadata = {}
        self.name = None
        self.radius_max = None
    
    def __str__(self):
        if self.data is not None and self.__nlayers__ > 0:
            output = "%s is a one-dimensional model with %s layers and radius up to %s km" % (self.name, self.__nlayers__,self.radius_max/1000.)
        else:
            output = "No model has been read into this reference1D instance yet"
        return output
            
    def read(self,file,fmt='card'):
        '''
        Read a card deck file used in OBANI. Other formats not ready yet
        '''
        if fmt=='card':
            names=['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta']
            formats=[np.float for ii in range(len(names))]
            modelarr = np.genfromtxt(file,dtype=None,comments='#',skip_header=3,
            names=names)
            self.metadata['attributes'] = names
            self.metadata['description'] = 'Read from '+file
            self.metadata['filename'] = file
            self.name = ntpath.basename(file)
        else:
            raise NotImplementedError('model format ',fmt,' is not currently implemented in reference1D.read')

        self.__nlayers__ = len(modelarr['radius'])
        # Create data array
        Model1D_Attr = np.dtype([(native_str(names[ii]),formats[ii]) for ii in range(len(names))])
        self.data = np.zeros(self.__nlayers__,dtype=Model1D_Attr)
        self.data['radius'] = modelarr['radius']
        self.data['rho'] = modelarr['rho']
        self.data['vpv'] = modelarr['vpv']
        self.data['vsv'] = modelarr['vsv']
        self.data['Qkappa'] = modelarr['Qkappa']
        self.data['Qmu'] = modelarr['Qmu']
        self.data['vph'] = modelarr['vph']
        self.data['vsh'] = modelarr['vsh']
        self.data['eta'] = modelarr['eta']
        self.radius_max = np.max(self.data['radius'])

    def get_Love_elastic(self):
        '''
        Get the Love parameters and Voigt averaged elastic properties with depth
        '''
        if self.data is not None and self.__nlayers__ > 0:
            self.data=append_fields(self.data, 'A', self.data['rho']*self.data['vph']**2 , usemask=False)
            self.data=append_fields(self.data, 'C', self.data['rho']*self.data['vpv']**2 , usemask=False)
            self.data=append_fields(self.data, 'N', self.data['rho']*self.data['vsh']**2 , usemask=False)
            self.data=append_fields(self.data, 'L', self.data['rho']*self.data['vsv']**2 , usemask=False)
            self.data=append_fields(self.data, 'F', self.data['eta']*(self.data['A']-2.*self.data['L']) , usemask=False)
            self.data=append_fields(self.data, 'kappa', (4.0*(self.data['A']+self.data['F']-self.data['N'])+self.data['C'])/9. , usemask=False)
            self.data=append_fields(self.data, 'mu', (self.data['A']+self.data['C']-2.*self.data['F']+5.*self.data['N']+6.*self.data['L'])/15. , usemask=False)
            self.data=append_fields(self.data, 'vp', np.sqrt(np.divide((self.data['kappa']+4.*self.data['mu']/3.),self.data['rho'])) , usemask=False)
            self.data=append_fields(self.data, 'vs', np.sqrt(np.divide(self.data['mu'],self.data['rho'])) , usemask=False)
            with np.errstate(divide='ignore', invalid='ignore'): # Ignore warning about dividing by zero
                xi = np.power(np.divide(self.data['vsh'],self.data['vsv']),2)
            self.data=append_fields(self.data, 'xi', xi , usemask=False)
            self.data=append_fields(self.data, 'phi', np.power(np.divide(self.data['vpv'],self.data['vph']),2) , usemask=False)
        else:
            raise ValueError('reference1D object is not allocated')

    def get_custom_parameter(self,parameters):
        '''
        Get the arrays of custom parameters defined in various Earth models
        '''
        if self.data is not None and self.__nlayers__ > 0:
            # convert to array for ease of looping
            if isinstance(parameters,string_types): parameters = np.array(parameters) 
            
            for parameter in parameters:
                if 'SH-SV' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['vsh'] - self.data['vsv'] , usemask=False)
                elif 'PH-PV' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['vph'] - self.data['vpv'] , usemask=False)
                elif '(SH+SV)*0.5' in parameter:
                    self.data=append_fields(self.data, parameter, (self.data['vsh'] + self.data['vsv'])/2. , usemask=False)
                elif '(PH+PV)*0.5' in parameter:
                    self.data=append_fields(self.data, parameter, (self.data['vph'] + self.data['vpv'])/2. , usemask=False)
                elif 'dETA/ETA' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['eta'] , usemask=False)
                elif 'dRHO/RHO' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['rho'] , usemask=False)                
                else:
                    raise NotImplementedError('parameter ',parameter,' is not currently implemented in reference1D.get_custom_parameter')
        else:
            raise ValueError('reference1D object is not allocated')

    def evaluate_at_depth(self,depth_in_km,parameter='vs',interpolation='linear'):   
        '''
        Get the values of a parameter at a given depth
        '''
        values=None
        if self.data is not None and self.__nlayers__ > 0:
            if parameter in self.data.dtype.names:
                values = self.data[parameter]
                depth_array = (6371000. - self.data['radius'])/1000. # in km
                # Sort to make interpolation possible
                indx = depth_array.argsort()
                # convert to array for ease of looping
                if isinstance(depth_in_km,float) or isinstance(depth_in_km,int): depth_in_km = np.array(depth_in_km) 
                values = griddata(depth_array[indx], values[indx], depth_in_km, method=interpolation)
            else:
                raise ValueError('parameter '+parameter+' not defined in array')
        else:
            raise ValueError('reference1D object is not allocated')
        return values

    def to_TauPmodel(self,fmt='tvel'):
        '''
        Writes a model file that is compatible with TauP.
        file format options 'tvel' and 'nd'.

        Note: TauP can't handle zero shear velocity in the ocean layer...
          To work around this, zero values an ocean layer will be written 
          as 1e-4.
        '''
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            f = open(model_name+'.tvel','w')
            f.write('{} - P\n'.format(model_name))
            f.write('{} - S\n'.format(model_name))

            for i in range(0,len(self.data)):
                f.write('{:2.4f}   {:2.4f}   {:2.4f}    {:2.4f}\n'.format(
                   (self.radius_max - self.data['radius'][::-1][i]) / 1000.0,
                   self.data['vp'][::-1][i] / 1000.0,
                   self.data['vs'][::-1][i] / 1000.0,
                   self.data['rho'][::-1][i] / 1000.0))       
            f.close()
        else:
            raise ValueError('reference1D object is not allocated')


    def to_axisem(self,anelastic=True,anisotropic=True):
        '''
         Write 1D model to be used as an external model in axisem
        '''
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            f = open(model_name+'.bm','w')
            n_discon = 0

            if anelastic:
                f.write('ANELASTIC     T\n')
            else:
                f.write('ANELASTIC     F\n')

            if anisotropic:
                f.write('ANISOTROPIC     T\n')
            else:
                f.write('ANISOTROPIC     F\n')

            f.write('UNITS      m\n')

            if anisotropic:
                f.write('COLUMNS   radius    rho    vpv    vsv    qka    qmu    vph    vsh    eta\n')

            for i in range(0,len(self.data)):
                f.write('{}    {}    {}    {}    {}    {}    {}    {}    {}\n'.format(
                self.data['radius'][::-1][i],
                self.data['rho'][::-1][i],
                self.data['vpv'][::-1][i],
                self.data['vsv'][::-1][i],
                self.data['Qkappa'][::-1][i],
                self.data['Qmu'][::-1][i],
                self.data['vph'][::-1][i],
                self.data['vsh'][::-1][i],
                self.data['eta'][::-1][i]) )

                if i < len(self.data)-1 and self.data['radius'][::-1][i] == self.data['radius'][::-1][i+1]:
                    depth_here = (self.radius_max - self.data['radius'][::-1][i]) / 1000.0 
                    n_discon += 1
                    f.write('#    Discontinuity {}, depth {:6.2f} km\n'.format(n_discon,depth_here))
        else:
            raise ValueError('reference1D object is not allocated')
            
    def plot(self,figuresize=[7,12],height_ratios=[2, 2, 1],ifshow=True,outfile='.png'):
        """ 
        Plot the cards array in a PREM like plot
        """
        depthkmarr = (6371000. - self.data['radius'])/1000. # in km
        #Set default fontsize for plots
        plots.updatefont(10)
        fig = plt.figure(1, figsize=(figuresize[0],figuresize[1]))
        gs = gridspec.GridSpec(3, 1, height_ratios=height_ratios) 
        fig.patch.set_facecolor('white')
        ax01=plt.subplot(gs[0])
        ax01.plot(depthkmarr,self.data['rho']/1000.,'k')
        ax01.plot(depthkmarr,self.data['vsv']/1000.,'b')
        ax01.plot(depthkmarr,self.data['vsh']/1000.,'b:')
        ax01.plot(depthkmarr,self.data['vpv']/1000.,'r')
        ax01.plot(depthkmarr,self.data['vph']/1000.,'r:')
        mantle=np.where( depthkmarr < 2891.)
        ax01.plot(depthkmarr[mantle],self.data['eta'][mantle],'g')
        ax01.set_xlim([0., 6371.])
        ax01.set_ylim([0, 14])
    
        majorLocator = MultipleLocator(2)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(1)
        ax01.yaxis.set_major_locator(majorLocator)
        ax01.yaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax01.yaxis.set_minor_locator(minorLocator)
        
        majorLocator = MultipleLocator(2000)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(1000)
        ax01.xaxis.set_major_locator(majorLocator)
        ax01.xaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax01.xaxis.set_minor_locator(minorLocator)
        ax01.set_ylabel('Velocity (km/sec), density (g/cm'+'$^3$'+') or '+'$\eta$')
        
        for para,color,xloc,yloc in [("$\eta$",'g',1500.,2.),("$V_S$",'b',1500.,7.8),("$V_P$",'r',1500.,13.5),("$\\rho$",'k',1500.,4.5),("$V_P$",'r',4000.,9.2),("$\\rho$",'k',4000.,12.5),("$V_S$",'b',5500.,4.5)]:
            ax01.annotate(para,color=color,
            xy=(3, 1), xycoords='data',
            xytext=(xloc/6371., yloc/14.), textcoords='axes fraction',
            horizontalalignment='left', verticalalignment='top')


        ax11=plt.subplot(gs[1])
        top1000km=np.where( depthkmarr < 1000.)
        ax11.plot(depthkmarr[top1000km],self.data['rho'][top1000km]/1000.,'k')
        ax11.plot(depthkmarr[top1000km],self.data['vsv'][top1000km]/1000.,'b')
        ax11.plot(depthkmarr[top1000km],self.data['vsh'][top1000km]/1000.,'b:')
        ax12 = ax11.twinx()
        ax12.plot(depthkmarr[top1000km],self.data['vpv'][top1000km]/1000.,'r')
        ax12.plot(depthkmarr[top1000km],self.data['vph'][top1000km]/1000.,'r:')
        ax11.plot(depthkmarr[top1000km],self.data['eta'][top1000km],'g')
        ax11.set_xlim([0., 1000.])
        ax11.set_ylim([0, 7])
        ax12.set_xlim([0., 1000.])
        ax12.set_ylim([-2, 12])        
        ax11.set_ylabel('Shear velocity (km/sec), density (g/cm'+'$^3$'+') or '+'$\eta$')
        ax12.set_ylabel('Compressional velocity (km/sec)')
        for para,color,xloc,yloc in [("$\eta$",'g',150.,1.),("$V_{S}$",'b',150.,4.3),("$V_{P}$",'r',120.,5.5),("$\\rho$",'k',150.,3.8)]:
            ax11.annotate(para,color=color,
            xy=(3, 1), xycoords='data',
            xytext=(xloc/1000., yloc/7.), textcoords='axes fraction',
            horizontalalignment='left', verticalalignment='top')
        ax12.set_yticks(np.arange(6, 14, step=2))
        majorLocator = MultipleLocator(200)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(100)
        ax11.xaxis.set_major_locator(majorLocator)
        ax11.xaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax11.xaxis.set_minor_locator(minorLocator)

        
        ax21=plt.subplot(gs[2], sharex=ax11)
        with np.errstate(divide='ignore', invalid='ignore'): # Ignore warning about dividing by zero
            anisoVs=(self.data['vsh']-self.data['vsv'])*200./(self.data['vsh']+self.data['vsv'])
        anisoVp=(self.data['vph']-self.data['vpv'])*200./(self.data['vph']+self.data['vpv'])
        ax21.plot(depthkmarr[top1000km],anisoVs[top1000km],'b')
        ax21.plot(depthkmarr[top1000km],anisoVp[top1000km],'r')
        ax21.set_ylim([0, 4])        
        majorLocator = MultipleLocator(1)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(0.5)
        ax21.yaxis.set_major_locator(majorLocator)
        ax21.yaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax21.yaxis.set_minor_locator(minorLocator)
        for para,color,xloc,yloc in [('Q'+'$_{\mu}$','k',400.,2.5),("$a_{S}$",'b',150.,3.7),("$a_{P}$",'r',100.,1.8)]:
            ax21.annotate(para,color=color,
            xy=(3, 1), xycoords='data',
            xytext=(xloc/1000., yloc/4.), textcoords='axes fraction',
            horizontalalignment='left', verticalalignment='top')


        ax22 = ax21.twinx()
        ax22.plot(depthkmarr[top1000km],self.data['Qmu'][top1000km],'k')
        ax21.set_xlabel('Depth (km)')
        ax21.set_ylabel("$V_P$"+' or '+"$V_S$"+' anisotropy (%)')
        ax22.set_ylabel('Shear attenuation Q'+'$_{\mu}$')
        ax22.set_ylim([0, 400])        
        majorLocator = MultipleLocator(100)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(50)
        ax22.yaxis.set_major_locator(majorLocator)
        ax22.yaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax22.yaxis.set_minor_locator(minorLocator)
        if ifshow: plt.show()
        plt.savefig(self.name+outfile)


#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys,os
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int,list,tuple

import numpy as np #for numerical analysis
import pdb    #for the debugger pdb.set_trace()
import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases

import struct
import xarray as xr
import traceback
import pandas as pd
####################### IMPORT REM3D LIBRARIES  #######################################
from .common import read3dmodelfile
#######################################################################################
class Realization(object):
    '''
    A class for the realization of a 3D model
    '''
    #########################       magic       ##########################
    def __init__(self,file=None):
        self.metadata = None
        self.data = None
        self._name = None
        self._type = None
        self._refmodel = None    
        self._infile = None    
        if file is not None:  success = self._read(file)

    def __str__(self):
        if self._name is not None:
            output = "%s is a three-dimensional model with %s as the reference model, read from %s of type %s" % (self._name,self._refmodel,self._infile,self._type)
        else:
            output = "No three-dimensional model has been read into this model3d instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._infile})'.format(self=self)

    #########################       decorators       ##########################
    @property
    def type(self):
        return self._type

    @property
    def name(self):
        return self._name

    @property
    def refmodel(self):
        return self._refmodel
    #########################       methods       #############################
        
    def _read(self,file):
        """
        Try reading the file into resolution/realization either as ascii, hdf5 or nc4
        """
        if (not os.path.isfile(file)): raise IOError("Filename ("+file+") does not exist")
        success = True
        try:# try ascii
            self._readascii(file)
        except:
            var1 = traceback.format_exc()
            try: # try nc4
                ds = xr.open_dataset(file)
                self._readnc4(ds)
                ds.close()
            except:
                try: #first close the dataset if opened with xarray above
                    ds.close()
                except NameError:
                    ds = None
                var2 = traceback.format_exc()
                print(var1)
                print(var2)
                success = False
        if success: self._infile = file

    def _readascii(self,modelfile):
        """
        Reads a standard 3D model file. maxkern is the maximum number of radial kernels
        and maxcoeff is the maximum number of corresponding lateral basis functions.
        resolution and realization are the indices for the resolution level
        and the realization from a model ensemble (usually 0 if a single file)
        """
        if (not os.path.isfile(modelfile)): raise IOError("Filename ("+modelfile+") does not exist")

        # read mean model
        model=read3dmodelfile(modelfile)
        # store relevant fields
        self._name = model['data']['name']
        self.data = model['data']['coef']
        self.metadata = model['metadata']
        self._type = 'ascii'
        self._refmodel = model['metadata']['refmodel']

    def _readnc4(self,ds):
        """
        Read netCDF4 file into a resolution and realization of model3D class.

        Input Parameters:
        -----------------

        ds: xarray Dataset handle

        """

        # Store in a dictionary
        metadata = {}
        for key in ds.attrs.keys(): metadata[key] = ds.attrs[key]
        metadata['nhorpar']=1
        metadata['null_model']=None
        metadata['shortcite']=ds.attrs['name']
        metadata['ityphpar']=np.array([3])
        metadata['typehpar']=np.array(['PIXELS'], dtype='<U40')
        # check the pixel size
        pixlat = np.unique(np.ediff1d(np.array(ds.latitude)))
        pixlon = np.unique(np.ediff1d(np.array(ds.longitude)))
        if not len(pixlat)==len(pixlon)==1: raise AssertionError('only one pixel size allowed in xarray')
        if not pixlat.item()==pixlon.item(): raise AssertionError('same pixel size in both lat and lon in xarray')
        metadata['hsplfile']=np.array([str(pixlat[0])+' X '+str(pixlat[0])], dtype='<U40')
        lenarr = len(ds.latitude)*len(ds.longitude)
        metadata['xsipix']=np.array([[pixlat[0] for ii in range(lenarr)]])
        metadata['xlapix'] = np.zeros([1,lenarr])
        metadata['xlopix'] = np.zeros([1,lenarr])
        # get data keys
        data_keys = []
        for key,_ in ds.data_vars.items(): data_keys.append(key)
        indx = 0
        idep = 0
        if len(ds.dims) == 3:
            _ , nlat, nlon = ds[data_keys[0]].shape
            for ilat in range(nlat):
                for ilon in range(nlon):
                    metadata['xlapix'][0][indx] = ds[data_keys[0]][idep,ilat,ilon].latitude
                    metadata['xlopix'][0][indx] = ds[data_keys[0]][idep,ilat,ilon].longitude
                    indx += 1
        else:
            raise ValueError('dimensions != 3')


        #depth differences and get depth extents
        depdiff = np.ediff1d(ds.depth)
        deptop = np.copy(ds.depth)
        depbottom = np.copy(ds.depth)

        for ii in range(len(depdiff)-2):
            deptop[ii] = deptop[ii] - (2.*depdiff[ii]-depdiff[ii+1])/2.
            depbottom[ii] = depbottom[ii] + (2.*depdiff[ii]-depdiff[ii+1])/2.
        for ii in range(len(depdiff),len(depdiff)-3,-1):
            deptop[ii] = deptop[ii] - (2.*depdiff[ii-1]-depdiff[ii-2])/2.
            depbottom[ii] = depbottom[ii] + (2.*depdiff[ii-1]-depdiff[ii-2])/2.

        ## create keys
        metadata['numvar']=len(data_keys)
        metadata['varstr']=np.array(data_keys, dtype='<U150')
        desckern = []; ivarkern = []; icount=0; coef=None
        for key in data_keys:
            icount = icount+1
            if 'topo' in key:
                depth_range = key.split('topo')[1]
                nlat, nlon = ds[key].shape
                descstring = u'{}, delta, {} km'.format(key,depth_range)
                try:
                    coef = np.vstack([coef,ds[key].values.flatten()])
                except:
                    coef = ds[key].values.flatten()
                desckern.append(descstring)
                ivarkern.append(icount)
            else:
                _ , nlat, nlon = ds[key].shape
                for ii,_ in enumerate(deptop):
                    descstring = u'{}, boxcar, {} - {} km'.format(key,deptop[ii],depbottom[ii])
                    try:
                        coef = np.vstack([coef,ds[key][ii].values.flatten()])
                    except:
                        coef = ds[key][ii].values.flatten()
                    desckern.append(descstring)
                    ivarkern.append(icount)
        metadata['desckern']=np.array(desckern, dtype='<U150')
        metadata['nmodkern']=len(desckern); metadata['ivarkern']=np.array(ivarkern)
        metadata['ihorpar']=np.ones(len(desckern),dtype = np.int)

        # get the number of coeffients
        metadata['ncoefhor']=np.array([lenarr])
        metadata['ncoefcum']=np.cumsum([metadata['ncoefhor'][ihor-1] for ihor in metadata['ihorpar']])

        # store to the object
        self._name = ds.attrs['name']
        self.data = pd.DataFrame(coef)
        self.metadata = metadata

        #rename the name field only if it is None
        self._type = 'netcdf4'
        self._refmodel = model['metadata']['refmodel']

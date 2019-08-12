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
from six import string_types # to check if variable is string using isinstance

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
        if file is not None:  success = self.read(file)

    def __str__(self):
        if self._name is not None:
            output = "%s is a three-dimensional model with %s as the reference model, read from %s of type %s" % (self._name,self._refmodel,self._infile,self._type)
        else:
            output = "No three-dimensional model has been read into this model3d instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._name})'.format(self=self)

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

    def read(self,file):
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

        # change None to string since it cannot be stored in netcdf
        for key in ds.attrs.keys(): metadata[key] = None if ds.attrs[key].lower() == 'none' else ds.attrs[key]
        for var in ds.data_vars:
            for key in ds[var].attrs.keys():
                val = ds[var].attrs[key]
                if isinstance(val,string_types):
                    if ds[var].attrs[key].lower() == 'none': ds[var].attrs[key] = None

        metadata['nhorpar']=1
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
        for key in ds.data_vars:
            if key != 'pixel_width': data_keys.append(key)
        indx = 0
        nlat = ds.dims['latitude']
        nlon = ds.dims['longitude']
        for ilat in range(nlat):
            metadata['xlapix'][0][indx:indx+nlon]=ds['latitude'].values[ilat]
            metadata['xlopix'][0][indx:indx+nlon]=ds['longitude'].values
            indx += nlon

        ## create keys
        metadata['numvar']=len(data_keys)
        metadata['varstr']=np.array(data_keys, dtype='<U150')

        # pre-allocate coef array. final shape is (n_depths, nlat*nlon).
        # n_depth for topo keys differ from others
        coef_ndepth=0
        coef_nlatnon=nlat*nlon
        for key in data_keys:
            if len(ds[key].dims) == 2:
                coef_ndepth=coef_ndepth+1
            elif len(ds[key].dims) == 3:
                ndepth , _, _ = ds[key].shape
                coef_ndepth=coef_ndepth+ndepth
            else:
                raise ValueError('dimension of key '+key+' cannot be anything except 2/3')
        coef=np.zeros([coef_ndepth,nlat*nlon])

        desckern = []; ivarkern = []; icount=0; idepth=0
        for key in data_keys:
            icount = icount+1
            if len(ds[key].dims) == 2:
                depth_range = ds[key].attrs['depth']
                nlat, nlon = ds[key].shape
                descstring = u'{}, delta, {} km'.format(key,depth_range)
                coef[idepth,:]=ds[key].values.flatten()
                idepth=idepth+1
                desckern.append(descstring)
                ivarkern.append(icount)
            else:
                ndepth , nlat, nlon = ds[key].shape
                for ii,deptop in enumerate(ds[key].attrs['start_depths']):
                    descstring = u'{}, boxcar, {} - {} km'.format(key,deptop,ds[key].attrs['end_depths'][ii])
                    desckern.append(descstring)
                    ivarkern.append(icount)
                coef[idepth:idepth+ndepth,:]=np.reshape(ds[key].values,[ndepth,nlat*nlon])
                idepth=idepth+ndepth

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
        self._refmodel = ds.attrs['refmodel']

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
import warnings
import struct
import xarray as xr
import traceback
import pandas as pd
from progressbar import progressbar

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
            self.readascii(file)
        except:
            var1 = traceback.format_exc()
            try: # try nc4
                ds = xr.open_dataset(file)
                self.readnc4(ds)
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

    def readascii(self,modelfile):
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

    def readnc4(self,ds):
        """
        Read netCDF4 file into a resolution and realization of model3D class.

        Input Parameters:
        -----------------

        ds: xarray Dataset handle

        """

        # Store in a dictionary
        metadata = {}
        metadata['attrs']={}

        # change None to string since it cannot be stored in netcdf, also stor attrbutes
        for key in ds.attrs.keys():
            if isinstance(ds.attrs[key],string_types):
                metadata[key] = None if ds.attrs[key].lower() == 'none' else ds.attrs[key]
            else:
                metadata[key] = ds.attrs[key]
        for var in ds.data_vars:
            for key in ds[var].attrs.keys():
                val = ds[var].attrs[key]
                if isinstance(val,string_types):
                    if ds[var].attrs[key].lower() == 'none': ds[var].attrs[key] = None
            metadata['attrs'][var] = ds[var].attrs

        metadata['nhorpar']=1
        metadata['shortcite']=ds.attrs['name']
        metadata['ityphpar']=np.array([3])
        metadata['typehpar']=np.array(['PIXELS'], dtype='<U40')
        # check the pixel size
        pixlat = np.unique(np.ediff1d(np.array(ds.latitude)))
        pixlon = np.unique(np.ediff1d(np.array(ds.longitude)))
        if len(pixlat)==len(pixlon)==1:
            if not pixlat.item()==pixlon.item(): warnings.warn('same pixel size in both lat and lon in xarray')
        else:
            warnings.warn('multiple pixel sizes have been found for xarray'+str(pixlat)+str(pixlon))

        metadata['hsplfile']=np.array([str(pixlat[0])+' X '+str(pixlat[0])], dtype='<U40')
        lenarr = len(ds.latitude)*len(ds.longitude)
        metadata['xsipix']= np.zeros([1,lenarr])
        metadata['xlapix'] = np.zeros([1,lenarr])
        metadata['xlopix'] = np.zeros([1,lenarr])
        # get data keys
        data_keys = []
        for key in ds.data_vars:
            if key != 'pixel_width': data_keys.append(key)
        indx = 0
        nlat = ds.dims['latitude']
        nlon = ds.dims['longitude']
        metadata['shape'] = nlat,nlon
        for ilat in range(nlat):
            metadata['xsipix'][0][indx:indx+nlon]=ds['pixel_width'].values[ilat]
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
                # the order needs to be in C so that it corresponds to ravel of a
                # (depth, lat,lon) array of coeffiecients to get modelarr
                # This also makes it consistent with how queries are made to KDtree in
                # buildtree3D

                coef[idepth:idepth+ndepth,:]=np.reshape(ds[key].values,[ndepth,nlat*nlon],order='C')
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

    def to_xarray(self,outfile=None,complevel=9, engine='netcdf4', writenc4 = False):
        '''
        write an xarrary dataset from a rem3d formatted ascii file

        Input parameters:
        ----------------

        ascii_file: path to rem3d format output file

        outfile: output netcdf file

        setup_file: setup file containing metadata for the model

        complevel, engine: options for compression in netcdf file

        writenc4: write a netcdf4 file, if True

        '''

        check = self.metadata['typehpar']=='PIXELS'
        if len(check) != 1 or not check[0]: raise IOError('cannot output netcdf4 file for horizontal parameterizations : ',self.metadata['typehpar'])
        if outfile == None and writenc4:
            if self._name.endswith('.nc4'):
                outfile = self._name
            else:
                outfile = self._name+'.'+self.metadata['kernel_set'].name+'.rem3d.nc4'

        # get sizes
        shape = self.metadata['shape']
        lat_temp = self.metadata['xlapix'][0]
        lon_temp = self.metadata['xlopix'][0]
        siz_temp = self.metadata['xsipix'][0]
        sortbylon = np.all(lon_temp == sorted(lon_temp))

        # unique values for reshaping
        if sortbylon:
            lon = np.unique(lon_temp)
            lat = np.unique(lat_temp)
            pxw = np.unique(siz_temp)
        else:
            arr=pd.DataFrame(np.vstack([lon_temp,lat_temp,siz_temp]).T,columns =['lon', 'lat', 'pxw'])
            arr = arr.sort_values(by=['lon','lat'])
            lon = pd.unique(arr['lon'])
            lat = pd.unique(arr['lat'])
            pxw = pd.unique(arr['pxw'])

        if not len(pxw)==1: raise warnings.warn('more than 1 pixel size in variable '+variable, pxw)

        # Get all depths
        kernel = self.metadata['kernel_set']
        depths = np.array([])
        for variable in self.metadata['varstr']:
            deptemp = kernel.pixeldepths(variable)
            if np.any(deptemp != None): depths = np.concatenate([depths,deptemp])
        alldepths= np.sort(np.unique(depths))

        # loop over variables
        data_vars={}

        # get the grid sizes stored
        pixel_array= np.reshape(arr['pxw'].values,
                            (len(lat),len(lon)),order='F')
        data_vars['pixel_width']=(('latitude', 'longitude'), pixel_array)

        # store every variable
        for variable in self.metadata['varstr']:
            deptemp = kernel.pixeldepths(variable)
            if np.any(deptemp == None): # 2D variable
                raise NotImplementedEddor('this feature is not available for 2D parameters '+variable)
            else:
                data_array = np.zeros((len(alldepths),len(lat),len(lon)))
                # update the kernel descriptions
                varind = np.where(self.metadata['varstr']==variable)[0][0]+1
                kerind = np.where(self.metadata['ivarkern']==varind)[0]
                if len(kerind) != len(deptemp): raise AssertionError('number of depths selected does not corrspong to the number of layers')
                # find depth indices
                dep_indx = np.searchsorted(alldepths,deptemp)

                values = self.data.iloc[kerind]
                for idep in progressbar(range(values.shape[0])):
                    if not sortbylon:
                        arr=pd.DataFrame(np.vstack([lon_temp,lat_temp,siz_temp,values.iloc[idep]]).T,columns =['lon', 'lat', 'pxw','values'])
                        arr = arr.sort_values(by=['lon','lat'])
                        valuerow = arr['values'].values
                    else:
                        valuerow = values.iloc[idep]
                    data_array[dep_indx[idep],:,:] = np.reshape(valuerow,
                                    (len(lat),len(lon)),order='F')
                data_vars[variable]=(('depth','latitude', 'longitude'), data_array)

        # get the grid sizes stored
        ds = xr.Dataset(data_vars = data_vars,
                        coords = {'depth':alldepths,'latitude':lat,'longitude':lon})

        # exclude some fields from dataframe attributes
        exclude = ['ityphpar','typehpar', 'hsplfile', 'xsipix', 'xlapix', 'xlopix','numvar', 'varstr', 'desckern', 'nmodkern', 'ivarkern','ihorpar', 'ncoefhor','ncoefcum', 'kernel_set','attrs']
        for field in self.metadata.keys():
            if field not in exclude:
                ds.attrs[field] = self.metadata[field]
        for variable in self.metadata['varstr']:
            for field in self.metadata['attrs'][variable].keys():
                ds[variable].attrs[field] = self.metadata['attrs'][variable][field]

        # write to file
        if outfile != None:
            print('... writing netcdf4 file .... ')
            # write to netcdf
            comp = {'zlib': True, 'complevel': complevel}
            encoding = {var: comp for var in ds.data_vars}
            # change None to string since it cannot be stored in
            for key in ds.attrs.keys():
                if ds.attrs[key] is None: ds.attrs[key] = 'None'
            for var in ds.data_vars:
                for key in ds[var].attrs.keys():
                    if ds[var].attrs[key] is None: ds[var].attrs[key] = 'None'
            ds.to_netcdf(outfile,engine=engine,encoding=encoding)
            print('... written netcdf4 file '+outfile)
        return ds

#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard AVNI format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys,os
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int,list,tuple

import numpy as np #for numerical analysis
import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases
from six import string_types # to check if variable is string using isinstance
import warnings
import struct
import xarray as xr
import traceback
import pandas as pd

####################### IMPORT AVNI LIBRARIES  #######################################
from .common import read3dmodelfile
from ..tools import convert_to_swp,convert2nparray
from .kernel_set import Kernel_set
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
        if file is not None:  self.read(file)

    def __str__(self):
        if self._name is not None:
            output = "%s is a three-dimensional model with %s as the reference model, read from %s of type %s" % (self._name,self._refmodel,self._infile,self._type)
        else:
            output = "No three-dimensional model has been read into this model3d instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._name})'.format(self=self)

    def __getitem__(self,key):
        """returns metadata from key"""
        return self.metadata[key]

    def __setitem__(self,key,data):
        """sets data to key"""
        self.metadata[key] = data


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

    @property
    def keys(self):
        return self.metadata.keys()

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
        # try to get the kernel set
        try:
            self['kernel_set'] = Kernel_set(self.metadata.copy())
            self.decode_mapping()
            self.decode_units()
            self.decode_scaling()
            self.scale_coefficients()
        except:
            warnings.warn('Warning: kernel_set could not be initialized in realization instance for '+file)
            pass
            #print(traceback.format_exc())

    def decode_units(self):
        if 'attrs' not in self.keys: self['attrs']={}
        kernel=self['kernel_set']

        # check if already available
        done=True
        if 'varstr' not in kernel.keys: raise ValueError('decoding units not possible without varstr initialized from model file and attributes.ini')
        for variable in kernel['varstr']: #loop over all variables, grabbing
            for type in ['unit','absolute_unit']:
                if variable in self['attrs'].keys():
                    done &= (type in self['attrs'][variable].keys())
                else:
                    done = False
        if done: return

        for key in ['unit','absolute_unit']:
            if key not in kernel.keys:
                raise ValueError('decoding units not possible without '+key+' initialized from model file and attributes.ini')
        for variable in kernel['varstr']: #loop over all variables, grabbing
            if variable not in self['attrs'].keys(): self['attrs'][variable]={}
            for type in ['unit','absolute_unit']:
                found = False
                units={}
                other_units=None
                for key in kernel[type].keys():
                    if key.replace(" ", "").lower() in variable.replace(" ", "").lower():
                        if found: raise ValueError('String '+key+' maps to multiple radial kernels. Please be more specific in attributes.ini')
                        found=True
                        self['attrs'][variable][type]=kernel[type][key]
                    if key.lower()=='other': other_units=kernel[type][key]
                if not found :
                    if other_units != None:
                        self['attrs'][variable][type]=other_units
                    else:
                        raise ValueError('no '+type+' decoded for variable '+variable+' based on attributes.ini')

        # Get units for mapped variables
        if 'mapping' in kernel.keys:
            for mapped in kernel['mapping']:
                searchstr = self['attrs'][mapped]['mapping']['variables'][0]
                for variable in kernel['varstr']: #loop over all variables
                    if searchstr in variable:
                        for type in ['unit','absolute_unit']:
                            self['attrs'][mapped][type]=self['attrs'][variable][type]

    def decode_symbols(self,searchstr,unique=False):
        searchstr = searchstr.split()
        symbols_all = ['-','+','X','x','*','^','/']
        variables = [x for x in searchstr if x not in set(symbols_all)]
        symbols = [x for x in searchstr if x not in set(variables)]

        scaling=[]
        expect=True
        for indx,val in enumerate(variables): # find scaling exclude last entry
            if expect:
                try:
                    scaling.append(float(val))
                    variables.pop(indx)
                    expect=False
                except ValueError:
                    if expect: scaling.append(1.)
                    expect=True

        # find unique entries corresponding to radial kernels
        if unique:
            temp=[]
            kernel=self['kernel_set']
            for val in variables: temp.append(kernel.search_radial(val,unique=True))
            variables = temp

        if not(len(variables)==len(scaling)==len(symbols)+1):
            raise ValueError("something wrong with decoding symbols, scaling and varianles from  "+searchstr)
        return variables,scaling,symbols

    def scale_coefficients(self):
        """
        scale model coefficients based on scaling
        """
        kernel = self['kernel_set']
        scaling=kernel.scaling
        radial_basis = kernel.data['radial_basis']
        for output in scaling:
            # find all kernels for this variable
            outfind = kernel.find_radial(output)
            basis_out = radial_basis[output]

            # checks
            input = scaling[output]['variables']
            if len(input) > 1: raise NotImplementedError('cannot scale more than one variable for now')
            infind = kernel.find_radial(input[0])
            basis_in = radial_basis[input[0]]

            scale_constant = scaling[output]['scaling']
            # now loop and check if radial basis is compatible
            for row,indxout in enumerate(outfind['index']):
                indxin = infind['index'][row]
                if basis_in[row] != basis_out[row]: raise ValueError('incompatible radial basis: '+basis_in[row]._name+' / '+basis_out[row]._name)

                # check that the lateral basis type is same
                ihorpar_in = self['ihorpar'][indxin]
                ihorpar_out = self['ihorpar'][indxout]
                if ihorpar_in != ihorpar_out:
                    warnings.warn('switching lateral parameterization of radial kernel, '+basis_out[row].name+' from '+self['typehpar'][ihorpar_out-1]+' to '+self['typehpar'][ihorpar_in-1])
                    self['ihorpar'][indxout]=ihorpar_in

                #scale the coefficients
                self.data.iloc[indxout] = self.data.iloc[indxin] * scale_constant

                # update cumulative number of coefficients
                self['kernel_set']['ncoefcum']=np.cumsum([self['ncoefhor'][ihor-1] for ihor in self['ihorpar']])

    def decode_scaling(self):
        kernel=self['kernel_set']
        if 'scaling' not in kernel.keys: return
        tempscaling={}
        if kernel['scaling'] is not None:
            for search in kernel['scaling']:
                variables,scaling,symbols = self.decode_symbols(kernel['scaling'][search],unique=True)

                # get the corresponding radial kernel
                radker = kernel.search_radial(search,unique=True)

                # store
                tempscaling[radker] = {'variables': variables,'scaling':scaling,'symbols':symbols}
        kernel['scaling']=tempscaling

    def decode_mapping(self):
        if 'attrs' not in self.keys: self['attrs']={}
        kernel=self['kernel_set']
        # decode units for mapped variables
        if 'mapping' not in kernel.keys: return
        for mapped in kernel['mapping']:
            # decode string
            variables,scaling,symbols = self.decode_symbols(kernel['mapping'][mapped])

            # Only = or - allowed for mapping
            for symb in np.unique(symbols):
                if symb not in ['-','+']: raise ValueError('Only + or - allowed for mapping')
            #store
            if mapped not in self['attrs'].keys(): self['attrs'][mapped]={}
            self['attrs'][mapped]['mapping']={'variables': variables,'constants':scaling,'symbols':symbols}

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
        write an xarrary dataset from a avni formatted ascii file

        Input parameters:
        ----------------

        ascii_file: path to avni format output file

        outfile: output netcdf file

        setup_file: setup file containing metadata for the model

        complevel, engine: options for compression in netcdf file

        writenc4: write a netcdf4 file, if True

        '''

        check = self['typehpar']=='PIXELS'
        if len(check) != 1 or not check[0]: raise IOError('cannot output netcdf4 file for horizontal parameterizations : ',self['typehpar'])
        if outfile == None and writenc4:
            if self._name.endswith('.nc4'):
                outfile = self._name
            else:
                outfile = self._name+'.'+self['kernel_set'].name+'.avni.nc4'

        # get sizes
        shape = self['shape']
        lat_temp = self['xlapix'][0]
        lon_temp = self['xlopix'][0]
        siz_temp = self['xsipix'][0]
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

        if not len(pxw)==1: warnings.warn('more than 1 pixel size in '+self._name)

        # Get all depths
        kernel = self['kernel_set']
        depths = np.array([])
        for variable in self['varstr']:
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
        for variable in self['varstr']:
            deptemp = kernel.pixeldepths(variable)
            if np.any(deptemp == None): # 2D variable
                # update the kernel descriptions
                varind = np.where(self['varstr']==variable)[0][0]+1
                kerind = np.where(self['ivarkern']==varind)[0]
                values = self.data.iloc[kerind]
                if not sortbylon:
                    arr=pd.DataFrame(np.vstack([lon_temp,lat_temp,siz_temp,values]).T,columns =['lon', 'lat', 'pxw','values'])
                    arr = arr.sort_values(by=['lon','lat'])
                    valuerow = arr['values'].values
                else:
                    valuerow = values.values
                data_array = np.reshape(valuerow,(len(lat),len(lon)),order='F')
                data_vars[variable]=(('latitude', 'longitude'), data_array)
            else:
                data_array = np.zeros((len(alldepths),len(lat),len(lon)))
                # update the kernel descriptions
                varind = np.where(self['varstr']==variable)[0][0]+1
                kerind = np.where(self['ivarkern']==varind)[0]
                if len(kerind) != len(deptemp): raise AssertionError('number of depths selected does not corrspong to the number of layers')
                # find depth indices
                dep_indx = np.searchsorted(alldepths,deptemp)

                values = self.data.iloc[kerind]
                for idep in range(values.shape[0]):
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
                ds.attrs[field] = self[field]
        exclude = ['sh'] #exclude spherical harmonic coefficients
        for variable in self['varstr']:
            for field in self['attrs'][variable].keys():
                if field not in exclude:
                    ds[variable].attrs[field] = self['attrs'][variable][field]

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

    def to_harmonics(self,lmax=40,variables=None):
        if variables == None: variables = self['varstr']
        variables = convert2nparray(variables)
        check = self['typehpar']=='PIXELS'
        if len(check) != 1 or not check[0]: raise IOError('cannot output harmonics for horizontal parameterizations : ',self['typehpar'])
        # Convert to numpy array
        namelist = ['latitude','longitude','value']
        formatlist = ['f8','f8','f8']
        dtype = dict(names = namelist, formats=formatlist)
        epixarr = np.zeros((self.data.shape[1]),dtype=dtype)
        epixarr['latitude'] = self['xlapix'][0]
        epixarr['longitude'] = self['xlopix'][0]

        coef=np.zeros((self['nmodkern'],(lmax+1)**2))
        for ivar,field in enumerate(variables):
            layers = np.where(self['ivarkern']==ivar+1)[0]
            print('... calculating spherical harmonic expansion of # '+str(ivar+1)+' / '+str(len(variables))+' : '+field+' , '+str(len(layers))+' layers')
            for ii,layer in enumerate(layers[:10]):
                epixarr['value'] = self.data.iloc[layer].to_numpy()
                shmatrix = convert_to_swp(epixarr,lmax=lmax)
                row=[]
                for ii in np.arange(len(shmatrix)):
                    l = shmatrix['l'][ii]; m = shmatrix['m'][ii]
                    if m==0:
                        row.append(shmatrix['cos'][ii])
                    else:
                        row.append(shmatrix['cos'][ii])
                        row.append(shmatrix['sin'][ii])
                coef[layer]=row
        self.harmonics = pd.DataFrame(coef)

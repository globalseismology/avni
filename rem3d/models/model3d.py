#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import numpy as np #for numerical analysis
import pdb    #for the debugger pdb.set_trace()
import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases
import matplotlib.pyplot as plt

from scipy import sparse
from configobj import ConfigObj
import re
from copy import deepcopy
import struct
import h5py
import xarray as xr
import traceback
import pandas as pd
if (sys.version_info[:2] < (3, 0)): input = raw_input
####################### IMPORT REM3D LIBRARIES  #######################################
from .. import tools   
from .. import plots 
from .. import constants
from .common import read3dmodelfile
from .kernel_set import kernel_set
#######################################################################################

# 3D model class
class model3d(object):
    '''
    A class for 3D reference Earth models used in tomography
    '''
    def __init__(self,file=None,**kwargs):
        self.metadata ={}
        self.data = {}
        self.name = None
        self.type = None
        self.refmodel = None
        self.description = None
        self.add_resolution(realization=True)
        if file is not None: 
            if (not os.path.isfile(file)): raise IOError("Filename ("+file+") does not exist")
            try:# try hdf5 for the whole ensemble
                success1 = True
                hf = h5py.File(file, 'r')
                if kwargs:
                    self.readhdf5(hf,**kwargs)
                else:
                    self.readhdf5(hf)
                self.description = "Read from "+file
                self.infile = file
                hf.close()
            except: # try netcdf or ascii for a single model
                try: #first close the hdf5 if opened with h5py above
                    hf.close()
                except NameError:
                    hf = None
                success1 = False
                var1 = traceback.format_exc()
                try:
                    if kwargs:
                        success2 = self.read(file,resolution=0,realization=0,**kwargs)
                    else:
                        success2 = self.read(file,resolution=0,realization=0)
                except:
                    print(var1)
            if not success1 and not success2: raise IOError('unable to read '+file+' as ascii, hdf5 or netcdf4')
            # try to get the kernel set
            for resolution in self.metadata.keys():
                try:
                    self.metadata[resolution]['kernel_set'] = kernel_set(self.metadata[resolution])
                except:
                    print('Warning: Kernelset could not initialized for'+str(resolution))
                        
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
        return realization

    def add_resolution(self,name=None,realization = False):
        """
        Added a resolution level to the object. num_realization is the number 
        of realizations of model coefficients within that object.
        """      
        num_resolution = len(self.data)
        self.metadata['resolution_'+str(num_resolution)] = {'name':name}
        self.data['resolution_'+str(num_resolution)] = {}
        if realization: 
            self.data['resolution_'+str(num_resolution)]['realization_0'] = {'name':None,'coef':None}
        return num_resolution

    def read(self,file,resolution=0,realization=0,**kwargs):
        """
        Try reading the file into resolution/realization either as ascii, hdf5 or nc4
        """
        if (not os.path.isfile(file)): raise IOError("Filename ("+file+") does not exist")
        success = True
        try:# try ascii
            if kwargs:
                self.readascii(file,resolution=resolution,realization=realization,**kwargs)
            else:
                self.readascii(file,resolution=resolution,realization=realization)
        except:
            var1 = traceback.format_exc()
            try: # try nc4
                ds = xr.open_dataset(file)
                if kwargs:
                    self.readnc4(ds,resolution=resolution,realization=realization,**kwargs)
                else:
                    self.readnc4(ds,resolution=resolution,realization=realization)   
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
        if success:
            self.description = "Read from "+file
            self.infile = file
        return success        
        

    def readascii(self,modelfile,resolution=0,realization=0,**kwargs):
        """
        Reads a standard 3D model file. maxkern is the maximum number of radial kernels
        and maxcoeff is the maximum number of corresponding lateral basis functions.
        resolution and realization are the indices for the resolution level
        and the realization from a model ensemble (usually 0 if a single file)
        """    
        if (not os.path.isfile(modelfile)): raise IOError("Filename ("+modelfile+") does not exist")
    
        # read mean model  
        if kwargs:
            model=read3dmodelfile(modelfile,**kwargs)
        else:
            model=read3dmodelfile(modelfile)
        
        self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['name'] = model['data']['name']
        self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['coef'] = model['data']['coef']

        self.metadata['resolution_'+str(resolution)] = model['metadata']
        #rename the name field only if it is None
        if self.name == None : self.name = model['data']['name']
        self.description = "Read from "+modelfile
        self.infile = modelfile
        self.type = 'rem3d'
        self.refmodel = model['metadata']['refmodel']
        
        return 
        
    def readnc4(self,ds,resolution=0,realization=0,**kwargs):
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
            ndepth, nlat, nlon = ds[data_keys[0]].shape
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
                ndepth, nlat, nlon = ds[key].shape
                for ii in range(len(deptop)):
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
        self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['name'] = ds.attrs['name']
        self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['coef'] = pd.DataFrame(coef)
        
        self.metadata['resolution_'+str(resolution)] = metadata
        
        #rename the name field only if it is None
        if self.name == None : self.name = metadata['name']
        self.type = 'rem3d'
        

    def readhdf5(self,hf,query=None):
        """
        Reads a standard 3D model file from a hdf5 file
        
        Input Parameters:
        ----------------
        
        query : if None, use the model available if only one is included. 
                Choose from query hf.keys() if multiple ones are available
                
        hf: hdf5 handle from h5py
        """    
        if query is None:
            if len(hf.keys()) == 1: 
                query = hf.keys()[0]
            else:
                str=''
                for str1 in hf.keys(): str=str+' , '+str1
                raise ValueError("... choose one from multiple query "+str)
        # read mean model   
        for name,value in hf[query].attrs.items(): 
            try:
                setattr(self, name,value) 
            except:
                setattr(self, name,None) 
                        
        # loop over resolution
        if len(self.data) < len(hf[query].keys()):
            # add more resolutions
            for ii in range(len(hf[query].keys()) - len(self.data)): self.add_resolution()
        for resolution in hf[query].keys():
            g1 = hf[query][resolution]
            if not g1.attrs['type']=='resolution': raise AssertionError()
            for name,value in g1.attrs.items(): 
                try:
                    self.metadata[resolution][name] = value
                except:
                    self.metadata[resolution][name] = None
            #now loop over every realization
            for case in g1.keys():
                g2 = hf[query][resolution][case]
                if not g2.attrs['type']=='realization': raise AssertionError()
                kerstr = self.metadata[resolution]['kerstr']
                key = self.name+'/'+kerstr+'/'+resolution+'/'+case
                self.data[resolution][case]['coef'] = pd.DataFrame(tools.io.load_numpy_hdf(hf,key))
                self.data[resolution][case]['name'] = g2.attrs['name']     
        return   
        
          
    def writehdf5(self, outfile = None, overwrite = False):
        """
        Writes the model object to hdf5 file
        """
        if outfile == None: outfile = self.infile+'.h5'
        if overwrite:
            hf = h5py.File(outfile, 'w')
        else:
            hf = h5py.File(outfile, 'a')
        g1 = hf.require_group(self.name)
        
        if self.name != None: g1.attrs['name']=self.name
        if self.type != None: g1.attrs['type']=self.type
        if self.refmodel != None: g1.attrs['refmodel']=self.refmodel
        if self.description != None: g1.attrs['description']=self.description
        
        for ires in range(len(self.data)):
            g2 = g1.require_group('resolution_'+str(ires))
            # write metadata for this resolution
            g2.attrs['type']='resolution'
            try:
                name = self.metadata['resolution_'+str(ires)]['name']
                g2.attrs['name']=name
            except:
                print('Warning: No name found for resolution '+str(ires))
            for key in self.metadata['resolution_'+str(ires)].keys():
                keyval = self.metadata['resolution_'+str(ires)][key]
                try:
                    if keyval.dtype.kind in {'U', 'S'}:
                        keyval = [n.encode('utf8') for n in keyval]
                        g2.attrs[key]=keyval
                except:
                    if keyval != None and key != 'kernel_set': 
                        g2.attrs[key]=keyval
                    
            #now loop over every realization
            for icase in range(len(self.data['resolution_'+str(ires)])):
                g3 = g2.require_group('realization_'+str(icase))
                g3.attrs['type']='realization'
                try:
                    name = self.data['resolution_'+str(ires)] ['realization_'+str(icase)]['name']
                    g3.attrs['name']=name
                except:
                    print('Warning: No name found for resolution '+str(ires)+', realization '+str(icase))
                # write the data array in the appropriate position
                #kerstr = self.metadata['resolution_'+str(ires)]['kerstr']
                key = self.name+'/resolution_'+str(ires)+'/realization_'+str(icase)
                out = tools.df2nparray(self.data['resolution_'+str(ires)] ['realization_'+str(icase)]['coef'])
                tools.io.store_numpy_hdf(hf,key,out)
                                
        hf.close()
        print('... written to '+outfile)

    def evaluate_at_point(self,latitude,longitude,depth_in_km,parameter='vs',resolution=0,realization=0,interpolated=False,tree=None,nearest=1): 
        """
        Evaluate the mode at a location (latitude, longitude,depth)
        
        Input Parameters:
        ----------------
        
        parameter: variable whose value will be returned
        
        interpolated: If True, use KDTree from a predefined grid. If False, evaluated 
                      exactly using kernel_set instance.
        """           
        if self.type != 'rem3d': raise NotImplementedError('model format ',self.type,' is not currently implemented in reference1D.coeff2modelarr')
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        
        # convert to numpy arrays
        latitude = tools.convert2nparray(latitude)
        longitude = tools.convert2nparray(longitude)
        depth_in_km = tools.convert2nparray(depth_in_km)
        resolution = tools.convert2nparray(resolution,int2float=False)
            
        #compute for each resolution
        for res in resolution:
            if not interpolated:
                # get the projection matrix
                project = self.calculateproj(latitude=latitude,longitude=longitude,depth_in_km=depth_in_km,parameter=parameter,resolution=res)
                modelarr = self.coeff2modelarr(resolution=res,realization=realization)
                predsparse = project['projarr']*modelarr
                values = predsparse.data
            else:
                if tree==None:
                    kerstr = self.metadata['resolution_'+str(res)]['kerstr']
                    treefile = kerstr+'.KDTree.3D.pkl'
                    #check that the horizontal param is pixel based
                    xlopix = self.metadata['resolution_'+str(res)]['xlopix'][0]
                    xlapix = self.metadata['resolution_'+str(res)]['xlapix'][0]
                    depths = self.getpixeldepths(res,parameter)
                    depth_in_km = np.array([])
                    for depth in depths: depth_in_km = np.append(depth_in_km,np.ones_like(xlopix)*depth)
                    xlapix = np.repeat(xlapix,len(depths))
                    xlapix = np.repeat(xlopix,len(depths))
                    tree = tools.tree3D(treefile,xlapix,xlapix,constants.R/1000. - depth_in_km)
                # get the interpolation
                values = tools.querytree3D(tree,latitude,longitude,depth_in_km,qpts_rad,vals,nearest=nearest)
                
        return values

    def getpixeldepths(self,resolution,parameter):
        typehpar = self.metadata['resolution_'+str(resolution)]['typehpar']
        if not len(typehpar) == 1: raise AssertionError('only one type of horizontal parameterization allowed')
        for type in typehpar: 
            if not type == 'PIXELS': raise AssertionError('for interpolation with tree3D')
        kernel_set = self.metadata['resolution_'+str(resolution)]['kernel_set']
        kernel_param = kernel_set.data['radial_basis'][parameter]
        depths = []
        for radker in kernel_param:
            depths.append((radker.metadata['depthtop']+radker.metadata['depthbottom'])/2.)
        return np.asarray(depths)
        
    def coeff2modelarr(self,resolution=0,realization=0):
        """
        Convert the coeffiecient matrix from the file to a sparse model array. Add
        up all the resolutions if a list is provided.
        
        realization : index of the set of coefficients in an ensemble. Default is 0
        as there is only one set of coefficients when read from a model file.
        
        resolution : list of resolutions to include the the modelarray
        """
        if self.type != 'rem3d': raise NotImplementedError('model format ',self.type,' is not currently implemented in reference1D.coeff2modelarr')
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        
        # convert to numpy arrays
        resolution = tools.convert2nparray(resolution,int2float=False)
            
        # Loop over resolution levels, adding coefficients 
        for ir in range(len(resolution)):
            try:
                coefficients =self.data['resolution_'+str(resolution[ir])]['realization_'+str(realization)]['coef']
                #if modelarr is already made use it
                try:
                    modelarr = self.data['resolution_'+str(resolution[ir])]['realization_'+str(realization)]['modelarr']
                except:
                    # Loop over all kernel basis
                    for ii in range(len(self.metadata['resolution_'+str(resolution[ir])]['ihorpar'])): 
                        # number of coefficients for this radial kernel
                        ncoef = self.metadata['resolution_'+str(resolution[ir])]['ncoefhor'][self.metadata['resolution_'+str(resolution[ir])]['ihorpar'][ii]-1]
                        # first radial kernel and first tesselation level
                        if ii == 0 and ir == 0: 
                            modelarr = coefficients.iloc[ii][:ncoef]                        
                        else:
                            modelarr = modelarr.append(coefficients.iloc[ii][:ncoef],ignore_index=True) 
                    # Convert to sparse matrix
                    # take the transpose to make dot products easy 
                    modelarr = sparse.csr_matrix(modelarr).transpose() 
                    self.data['resolution_'+str(resolution[ir])]['realization_'+str(realization)]['modelarr'] = modelarr 
            except AttributeError: 
                raise ValueError('resolution '+str(resolution[ir])+' and realization '+str(realization)+' not filled up yet.')
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
            
    def readprojbinary(self,lateral_basis):
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
            tools.io.store_numpy_hdf(h5f,'grid',grid)
            h5f.attrs["model"] = np.string_(model)
            h5f.attrs["param"] = np.string_(lateral_basis)
            
            for ii in np.arange(len(deptharr)):
                for jj in np.arange(len(refstrarr)):
                    proj = h5f.create_group("projection/depth_"+str(ii)+ "/refstr_"+str(jj))
                    proj.attrs['refvalue']=refvalarr[ii,jj]
                    tools.io.store_sparse_hdf(h5f,"projection/depth_"+str(ii)+ "/refstr_"+str(jj),projarr[ii,jj])
            h5f.close()   
        else:
            print("....reading "+outfile)
            h5f = h5py.File(outfile,'r')
            ndp=h5f.attrs['ndp']; npx=h5f.attrs['npx']; nhorcum=h5f.attrs['nhorcum']
            neval=h5f.attrs['neval']; deptharr=h5f.attrs['deptharr']
            refstrarr=h5f.attrs['refstrarr']; grid = tools.io.load_numpy_hdf(h5f,'grid') 
            xlat=grid['xlat']; xlon=grid['xlon']; area=grid['area']
            model=h5f.attrs['model']; param=h5f.attrs['param']
            
            # Always read the following from hdf5 so projarr is list ofsparsehdf5-fast I/O
            projarr={};refvalarr={}
            for ii in np.arange(len(deptharr)):
                for jj in np.arange(len(refstrarr)):
                    proj = h5f["projection/depth_"+str(ii)+ "/refstr_"+str(jj)]
                    refvalarr[ii,jj] = proj.attrs['refvalue']
                    projarr[ii,jj] = tools.io.load_sparse_hdf(h5f,"projection/depth_"+str(ii)+ "/refstr_"+str(jj))
            h5f.close()  
     
        # Write to a dictionary
        projection = {}
        projection['ndp']=ndp; projection['npx']=npx; projection['nhorcum']=nhorcum; 
        projection['neval']=neval; projection['deptharr']=deptharr; projection['refstrarr']=refstrarr
        projection['xlat']=xlat; projection['xlon']=xlon; projection['area']=area
        projection['refvalarr']=refvalarr; projection['projarr']=projarr       
        projection['model']=model; projection['param']=lateral_basis         
        return projection

    def calculateproj(self,latitude,longitude,depth_in_km,parameter='(SH+SV)*0.5',resolution=0):
        """
        Get the projection matrix from a lateral basis to another and for particular depths  
        """    
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        
        # convert to numpy arrays
        latitude = tools.convert2nparray(latitude)
        longitude = tools.convert2nparray(longitude)
        depth_in_km = tools.convert2nparray(depth_in_km)
        parameter = tools.convert2nparray(parameter)
                   
        if not len(latitude)==len(longitude)==len(depth_in_km): raise AssertionError('latitude, longitude and depth_in_km should be of same length')
                               
        # Get the radial projection file
        kernel = self.metadata['resolution_'+str(resolution)]['kernel_set']

        # loop through parameter and append the projection for each location
        refstrarr=[]
        for param in parameter:
            for iloc in range(len(latitude)):
                lat = latitude[iloc]
                lon = longitude[iloc]
                dep = depth_in_km[iloc]
                if iloc == 0:
                    projarr = kernel.getprojection(lat,lon,dep,param)
                else:
                    projarr = sparse.vstack([projarr,kernel.getprojection(lat,lon,dep,param)])
                refstrarr.append(param)
        
        # Write to a dictionary
        projection = {}
        projection['projarr']=projarr 
        projection['ndp']=len(depth_in_km)
        projection['deptharr']=depth_in_km
        projection['refstrarr']=np.array(refstrarr)
#         projection['ndp']=ndp; projection['npx']=npx; projection['nhorcum']=nhorcum; 
#         projection['neval']=neval; projection['deptharr']=deptharr; projection['refstrarr']=refstrarr
#         projection['xlat']=xlat; projection['xlon']=xlon; projection['area']=area
#         projection['refvalarr']=refvalarr;       
#         projection['model']=model; projection['param']=lateral_basis         
        return projection


    def projslices(self,projection,variable,depth,resolution=0,realization=0):
        """
        Projection matrix multiplied by model ensemble. Choses the nearest depth available for the projection.
        """    
        if self.name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")

        modelarr = self.coeff2modelarr(resolution=resolution,realization=realization)
        # Select appropriate arrays from projection matrix, read from file
        projarr = projection['projarr']; refstrarr = projection['refstrarr']; deptharr = projection['deptharr']  
    
        varindex = np.where(refstrarr==variable)[0]
        if len(varindex) == 0: raise ValueError('no projection found for variable '+variable)
        absdiff=abs(deptharr-depth)
        # Find the nearest depth
        depindex = absdiff.argmin()
        if min(absdiff) > 0. :
            print ("No unique depth found in the projection matrix. Choosing the nearest available depth "+str(deptharr[depindex]))
        modelselect=projarr[depindex,varindex[0]]*modelarr    
        return modelselect,deptharr[depindex]
        
    def reparameterize(self, model3d,resolution=0,realization=0):
        """
        Inverts for new coefficients in self.data from the coefficients in model3d class
        """
        if type(model3d).__name__ != 'model3d': raise ValueError('input model class should be a model3d instance')
        
        # Get the radial projection file
        selfmeta = self.metadata['resolution_'+str(resolution)]
        newmeta = model3d.metadata['resolution_'+str(resolution)]
        selfkernel = selfmeta['kernel_set']
        newkernel = newmeta['kernel_set']
        
        # get the projection matrix for each variable in self
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        for variable in selfmeta['varstr']:
            ivarfind =np.where(selfmeta['varstr']==variable)[0]
            if not len(ivarfind) == 1: raise AssertionError('only one parameter can be selected in eval_kernel_set')
            findvar = selfmeta['varstr'][ivarfind[0]]
            findrad = np.array([(ii, selfmeta['desckern'][ii]) for ii in np.arange(len(selfmeta['ivarkern'])) if ivarfind[0]+1 == selfmeta['ivarkern'][ii]],dtype=dt)
            
            # find the corresponding radial kernels in newmeta
            ivarfind2 =np.where(newmeta['varstr']==variable)[0]
            pdb.set_trace()
            if len(ivarfind2) != 1: 
                stringout = ''
                for ii in range(len(newmeta['varstr'])): stringout=stringout+str(ii)+'. '+newmeta['varstr'][ii]+'\n'
                print('')
                print(stringout)
                try:
                    x = int(input('Warning: no unique corresponding variable found for '+variable+'. Select one index from above: - default is 1:'))
                except (ValueError,EOFError):
                    x = 1
                ivarfind2 = np.array([x])
            findvar2 = newmeta['varstr'][ivarfind2[0]]
            findrad2 = np.array([(ii, newmeta['desckern'][ii]) for ii in np.arange(len(newmeta['ivarkern'])) if ivarfind2[0]+1 == newmeta['ivarkern'][ii]],dtype=dt)
            
            pdb.set_trace()
            # calculate the projection matrix for all locations
            longitude = selfmeta['xlopix'][0]; latitude = selfmeta['xlapix'][0]
            radialinfo = selfkernel.data['radial_basis'][variable][0].metadata
            depth_in_km = np.average(np.array([radialinfo['depthtop'],radialinfo['depthbottom']]),axis=0)
            # loop over depths and append the projection matrices
            proj = model3d.calculateproj(latitude=np.tile(latitude,len(depth_in_km)),longitude=np.tile(longitude,len(depth_in_km)),depth_in_km=np.repeat(depth_in_km,len(latitude)),parameter=findvar2,resolution=resolution)
                
            pdb.set_trace()
       

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

    def plotslices(self,lateral_basis='pixel1',dbs_path='~/dbs',x=0,percent_or_km='%',colormin = -6.,colormax=6.,depth=None,resolution=0,realization=0):
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
        if not isinstance(resolution, int): raise TypeError('resolution must be an integer, not %s' % type(resolution))
        if not isinstance(realization, int): raise TypeError('realization must be an integer, not %s' % type(realization))

        
        typehpar = self.metadata['resolution_'+str(resolution)]['typehpar']
        if len(typehpar) != 1 or typehpar[0] != 'PIXELS': raise ValueError('Slices can only be made for pixel paramterization')
            
        # Select appropriate arrays from projection matrix, read from file        
        lat = self.metadata['resolution_'+str(resolution)]['xlapix'][0]
        lon = self.metadata['resolution_'+str(resolution)]['xlopix'][0]

        refstrarr = self.metadata['resolution_'+str(resolution)]['varstr']
        # select models based on parameter and depth desired
        new_figure='y'  # flag for done
        colormin = -6.
        colormax = 6.
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
                ifplot =True
                try:
                    for ii in np.arange(len(refstrarr)): print(ii,refstrarr[ii])
                    try:
                        x = int(input("Select variable to plot - default is 0:"))
                    except (ValueError,EOFError):
                        x = x
                    
                    if 'topo' in refstrarr[x]:
                        #find the radial kernels for this paramter
                        kerfind = np.where(self.metadata['resolution_'+str(resolution)]['ivarkern']==x+1)[0]
                        if len(kerfind) == 1:
                            modelarray = self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['coef'].iloc[kerfind[0]]
                        else:
                            flag=flag-1
                            ifplot =False
                    else:
                        # get the depths available for this parameter
                        deptharr = self.getpixeldepths(resolution,refstrarr[x])
                        #depth differences and get depth extents
                        depdiff = np.ediff1d(deptharr)
                        deptop = np.copy(deptharr)
                        depbottom = np.copy(deptharr)
                        for ii in range(len(depdiff)-2):
                            deptop[ii] = deptop[ii] - (2.*depdiff[ii]-depdiff[ii+1])/2.
                            depbottom[ii] = depbottom[ii] + (2.*depdiff[ii]-depdiff[ii+1])/2.
                        for ii in range(len(depdiff),len(depdiff)-3,-1):
                            deptop[ii] = deptop[ii] - (2.*depdiff[ii-1]-depdiff[ii-2])/2.
                            depbottom[ii] = depbottom[ii]+ (2.*depdiff[ii-1]-depdiff[ii-2])/2.
                    
                        try:
                            depth = float(input("Select depth - select any value for topography ["+str(round(min(deptop),2))+"-"+str(round(max(depbottom),2))+"] :"))
                        except (ValueError,EOFError):
                            if depth is None:
                                depth = min(deptharr)
                            else:
                                depth = depth
                        if depth < min(deptop) or depth > max(depbottom):
                            flag=flag-1
                            ifplot =False
                        else:
                            #find the radial kernels for this paramter
                            kerfind = np.where(self.metadata['resolution_'+str(resolution)]['ivarkern']==x+1)[0]
                            #evaluate at all points
                            ind = np.where(np.logical_and(depth>deptop, depth<=depbottom))[0][0]
                            modelarray = self.data['resolution_'+str(resolution)]['realization_'+str(realization)]['coef'].iloc[kerfind[ind]]

                    if ifplot:
                        # Get limits for colorbar
                        try: 
                            colorstr = input("Input two values for minimum and maximum values of colorbar - default is "+str(colormin)+" "+str(colormax)+":")
                            colormin,colormax = float(colorstr.split()[0]),float(colorstr.split()[1])
                        except (ValueError,IndexError,EOFError):
                            colormin = colormin; colormax=colormax
            
                        # Plot the model
                        test = np.vstack((lat,lon,modelarray)).transpose()
                        dt = {'names':['lat', 'lon', 'val'], 'formats':[np.float, np.float, np.float]}
                        plotmodel = np.zeros(len(test), dtype=dt)
                        plotmodel['lat'] = test[:,0]; plotmodel['lon'] = test[:,1]; plotmodel['val'] = test[:,2]
                        ax=fig.add_subplot(subploty,subplotx,flag)
                        plots.globalmap(ax,plotmodel,colormin,colormax,dbs_path=dbs_path, colorlabel='Anomaly', grid=[30.,90.],gridwidth=0,projection='robin',lat_0=0, lon_0=150., colorpalette='rem3d',colorcontour=21)
                        ax.set_title(refstrarr[x]+' at '+str(depth)+' km.' if 'Topo' not in refstrarr[x] and 'topo' not in refstrarr[x] else refstrarr[x])
                        fig.canvas.draw() 
                except SyntaxError:
                    flag=flag-1
            try:
                new_figure = input("Another figure? y/n:")
            except (EOFError):
                new_figure = 'n'
    
        return 


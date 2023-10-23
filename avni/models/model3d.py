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
import pdb    #for the debugger pdb.set_trace()
import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases

from scipy import sparse,linalg
import re
import struct
import h5py
import traceback
import pandas as pd
import copy
import warnings
from six import string_types # to check if variable is string using isinstance
import pint
import xarray as xr

####################### IMPORT AVNI LIBRARIES  #######################################
from .. import tools
from .. import constants
from ..mapping import spher2cart,inpolygon
from .kernel_set import Kernel_set
from .realization import Realization
from .reference1d import Reference1D
from .common import getLU2symmetric,readResCov
from .profiles import Profiles
from ..f2py import drspln,drsple # cubic spline interpolation used in our codes
from ..plots import plot1section

#######################################################################################

# 3D model class
class Model3D(object):
    '''
    A class for 3D reference Earth models used in tomography
    '''
    #########################       magic       ##########################
    def __init__(self,file=None,**kwargs):
        self.metadata ={}
        self.data = {}
        self._name = None
        self._type = None
        self._infile = None
        self._refmodel = None
        self._description = None
        if file is not None: self.read(file,**kwargs)

    def __str__(self):
        if self._name is not None:
            output = "%s is a three-dimensional model ensemble with %s resolution levels, %s realizations of %s type and %s as the reference model" % (self._name,self.num_resolutions, len(self.data['resolution_0']),self._type,self._refmodel)
        else:
            output = "No three-dimensional model has been read into this model3d instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._name})'.format(self=self)

    def __len__(self):
        return len(self.metadata)

    def __getitem__(self,key):
        """returns data and metadata from key, a tuple of (resolution, realization)"""
        key = tools.convert2nparray(key,int2float=False)
        # if only first index, return metadata for the resolution,
        if not len(key) in [1,2] : raise AssertionError()
        metadata = self.metadata['resolution_'+str(key[0])]
        if len(key) == 1:
            return metadata
        else:
            data = self.data['resolution_'+str(key[0])]['realization_'+str(key[1])]
            return data

    def __setitem__(self,key,data):
        """sets (data, metadata) data_meta to key, a tuple of (resolution, realization)"""
        if isinstance(key, (int,np.int64)):
            self.metadata['resolution_'+str(key)] = data
        elif isinstance(key, tuple):
            self.data['resolution_'+str(key[0])]['realization_'+str(key[1])] = data
        else:
            raise TypeError('invalid input type for model3d')

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
            setattr(result, k, copy.deepcopy(v, memo))
        return result

    #########################       decorators       ##########################

    @property
    def name(self):
        return self._name

    @property
    def num_resolutions(self):
        return len(self.metadata)

    @property
    def num_realizations(self):
        output = []
        for res in self.metadata: output.append(len(self.data[res]))
        return output

    #########################       methods       #############################

    def read(self,file,**kwargs):
        if (not os.path.isfile(file)): raise IOError("Filename ("+file+") does not exist")
        success = False
        try:# try hdf5 for the whole ensemble
            hf = h5py.File(file, 'r')
            if kwargs:
                self.readhdf5(hf,**kwargs)
            else:
                self.readhdf5(hf)
            self._description = "Read from "+file
            self._infile = file
            hf.close()
            success = True
        except: # try netcdf or ascii for a single model
            try: #first close the hdf5 if opened with h5py above
                hf.close()
            except NameError:
                hf = None
            var1 = traceback.format_exc()
            try:
                realization = Realization(file)
                # store in the last resolution
                self.add_resolution(metadata=realization.metadata)
                self.add_realization(coef=realization.data,name=realization._name)
                self._name = realization._name
                self._type = realization._type
                self._refmodel = realization._refmodel
                self._infile = file
                success = True
            except:
                print('############    Tried reading as hdf5   ############')
                print(var1)
                print('############    Tried reading as ascii   ############')
                print(traceback.format_exc())
        if not success: raise IOError('unable to read '+file+' as ascii, hdf5 or netcdf4')

    def plot(self):
        from ..plots import plotmodel3d
        plotmodel3d(self)

    def add_realization(self,coef=None,name=None,resolution=None):
        """
        Added a set of realizations to the object. resolution is the tesselation level at
        which the instance is update with name and coeff (default: None, last resolution).
        """
        if resolution==None:
            if self.num_resolutions == 0: raise ValueError('add_resolution first') # add a resolution if none
            resolution = self.num_resolutions - 1 # last available resolution
            realization = self.num_realizations[resolution]
        else:
            if isinstance(resolution, int) :
                realization = self.num_realizations['resolution_'+str(resolution)]
            else:
                raise ValueError('Invalid resolution in add_realization')
        #initialize
        self.data['resolution_'+str(resolution)]['realization_'+str(resolution)] = {}
        # fill it up
        if name == None: name = str(realization)
        data = {'name':name,'coef':coef}
        self[resolution,realization]=data
        return

    def add_resolution(self,metadata=None):
        """
        Added a resolution level to the object. num_realizations is the number
        of realizations of model coefficients within that object.
        """
        current = self.num_resolutions
        if metadata == None: metadata = {'name':str(current)}
        self[current] = metadata
        self.data['resolution_'+str(current)] = {} # empty resolutions
        return

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
                strsum=''
                for str1 in hf.keys(): strsum=strsum+' , '+str1
                raise ValueError("... choose one from multiple query "+strsum)
        # read mean model
        for name,value in hf[query].attrs.items():
            try:
                setattr(self, name,value)
            except:
                setattr(self, name,None)

        # loop over resolution
        if len(self.data) < len(hf[query].keys()):
            # add more resolutions
            for _ in range(len(hf[query].keys()) - len(self.data)): self.add_resolution()
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
                key = self._name+'/'+kerstr+'/'+resolution+'/'+case
                self.data[resolution][case]['coef'] = pd.DataFrame(tools.io.load_numpy_hdf(hf,key))
                self.data[resolution][case]['name'] = g2.attrs['name']
        return

    def writerealization(self, resolution=0,realization=0,outfile=None):
        realize = Realization()
        # transfer data and metadata
        realize.metadata = self[resolution]
        realize._name = self._name
        realize.data = self[resolution,realization]['coef']
        realize._type = self._type
        realize._refmodel = self._refmodel
        realize._infile = self._infile

        #write the file
        realize.to_xarray(outfile=outfile,writenc4=True)

    def writehdf5(self, outfile = None, overwrite = False):
        """
        Writes the model object to hdf5 file
        """
        if outfile == None: outfile = self._infile+'.h5'
        if overwrite:
            hf = h5py.File(outfile, 'w')
        else:
            hf = h5py.File(outfile, 'a')
        g1 = hf.require_group(self._name)

        if self._name != None: g1.attrs['name']=self._name
        if self._type != None: g1.attrs['type']=self._type
        if self._refmodel != None: g1.attrs['refmodel']=self._refmodel
        if self._description != None: g1.attrs['description']=self._description

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
                key = self._name+'/resolution_'+str(ires)+'/realization_'+str(icase)
                out = tools.df2nparray(self.data['resolution_'+str(ires)] ['realization_'+str(icase)]['coef'])
                tools.io.store_numpy_hdf(hf,key,out)

        hf.close()
        print('... written to '+outfile)

    def buildtree3D(self,resolution=0,dbs_path=tools.get_filedir()):
        """
        Build a KDtree interpolant based on the metadata

        """
        typehpar = self[resolution]['typehpar']
        if not len(typehpar) == 1: raise AssertionError('only one type of horizontal parameterization allowed')
        for types in typehpar:
            if not types == 'PIXELS': raise AssertionError('for interpolation with tree3D')

        # not make list of depths
        kerstr = self[resolution]['kernel_set'].name
        check = True; depth_shared=None
        for indx,variable in enumerate(self[resolution]['varstr']):
            depths = self.getpixeldepths(resolution,variable)
            if not np.any(depths == None): # ignore 2D variables in this comparisons
                if np.any(depth_shared == None):
                    depth_shared = depths
                else:
                    check = check and (np.all(depth_shared==depths))
        if not check: raise ValueError('All 3D variables need to share the same radial parameterization for a common KDtree to work')
        #get full path
        dbs_path=tools.get_fullpath(dbs_path)
        treefile = dbs_path+'/'+constants.planetpreferred+'.'+kerstr+'.KDTree.3D.pkl'
        if os.path.isfile(treefile):
            tree = tools.tree3D(treefile)
        else:
            #check that the horizontal param is pixel based
            xlapix = self[resolution]['xlapix'][0]
            xlopix = self[resolution]['xlopix'][0]
            nlat = len(xlapix)
            ndep = len(depth_shared)
            depth_in_km = np.zeros(nlat*ndep)
            for indx,depth in enumerate(depth_shared):
                depth_in_km[indx*nlat:(indx+1)*nlat] = depth * np.ones(nlat)
            xlapix = np.tile(xlapix,ndep)
            xlopix = np.tile(xlopix,ndep)
            tree = tools.tree3D(treefile,xlapix,xlopix,constants.R.to('km').magnitude - depth_in_km)
        return tree


    def ifwithinregion(self,latitude,longitude,depth_in_km=None,resolution=0):
        # check if the queries is within the bounds of the model
        lat_max = float(self[resolution]['geospatial_lat_max'])
        lat_min = float(self[resolution]['geospatial_lat_min'])
        lon_min = float(self[resolution]['geospatial_lon_min'])
        lon_max = float(self[resolution]['geospatial_lon_max'])
        lon_min = lon_min + 360. if lon_min < 0. else lon_min
        lon_max = lon_max + 360. if lon_max < 0. else lon_max
        dep_min = self[resolution]['geospatial_vertical_min']
        dep_max = self[resolution]['geospatial_vertical_max']
        # convert to numpy arrays
        latitude = tools.convert2nparray(latitude)
        longitude = tools.convert2nparray(longitude)
        if np.any(depth_in_km == None):
            if not(len(latitude)==len(latitude)): raise AssertionError('latitude and longitude need to be of same size')
            checks = (latitude <= lat_max) & (latitude >= lat_min) & \
                 (longitude <= lon_max) & (longitude >= lon_min) & \
                 (depth_in_km <= dep_max) & (depth_in_km >= dep_min)
        else:
            depth_in_km = tools.convert2nparray(depth_in_km)
            if not(len(latitude)==len(latitude)==len(depth_in_km)): raise AssertionError('latitude, longitude and depth need to be of same size')
            checks = (latitude <= lat_max) & (latitude >= lat_min) & \
                 (longitude <= lon_max) & (longitude >= lon_min)
        return checks

    def average_in_polygon(self,depth_in_km,parameter,polygon_latitude, polygon_longitude,num_cores=1, orientation = 'anti-clockwise', threshold = 1E-6,grid=1.,outside=False,mask=None,**kwargs):
        """
        Get the average values within a polygon

        Input Parameters:
        ----------------

        parameter: variable whose value will be returned

        depth_in_km: depth where values are being queried at (km)

        polygon_latitude,polygon_longitude: closed points that define the polygon.
                                First and last points need to be the same.

        grid: fineness of the grid to use for averaging (in degrees)

        num_cores: Number of cores to use for the calculations

        orientation: clockwise or anti-clockwise orientation of points specified above

        threshold: limit to which the sum of azimuth check to (-)360 degrees is permitted
               to be defined as within the polygon.

        outside: give average values outside polygon

        Output:
        ------

        average: average value within the polygon of interest

        tree: if kwarg argument to evaluate_at_location has interpolated=True and
                tree==None, returns the tree data structure.
        """

        ## convert to numpy arrays. Various checks
        polylat = tools.convert2nparray(polygon_latitude)
        polylon = tools.convert2nparray(polygon_longitude)
        if len(polylat) != len(polylon): raise ValueError('polygon latitude and longitude need to be of same length')
        if (polygon_latitude[-1] != polygon_latitude[0]) or (polygon_longitude[-1] != polygon_longitude[0]): raise ValueError('polygon should be closed and therefore have same start and end points')
        nvertices = len(polylon)
        if orientation not in ['clockwise','anti-clockwise']: raise ValueError('orientation needs to be clockwise or anti-clockwise')

        ## define the lat/lon options based on grid
        if isinstance(grid,(int,np.int64,float)):
            latitude = np.arange(-90+grid/2., 90,grid)
            longitude = np.arange(0+grid/2., 360,grid)
            meshgrid =  True
            nlat = len(latitude)
            nlon = len(longitude)
            nprofiles = nlat*nlon
        elif isinstance(grid,string_types):
            meshgrid =  False
            raise NotImplementedError('needs to be implemented soon')

        ## Define the output array and area
        data_vars={}
        data_vars['inside']=(('latitude', 'longitude'), np.zeros((nlat,nlon)))
        data_vars['mask']=(('latitude', 'longitude'), np.zeros((nlat,nlon), dtype=bool))
        if outside:data_vars['outside']=(('latitude', 'longitude'), np.zeros((nlat,nlon)))
        outarr = xr.Dataset(data_vars = data_vars,
                coords = {'latitude':latitude,'longitude':longitude})
        area = tools.areaxarray(outarr)

        ## Loop over all items, checking if it is inside the column
        if mask is not None:
            outarr['mask'] = mask
        else:
            for ii,lat in enumerate(outarr.coords['latitude'].values):
                for jj,lon in enumerate(outarr.coords['longitude'].values):
                    within = inpolygon(lat,lon,polylat,polylon,num_cores,orientation,threshold)
                    if within: outarr['mask'][ii,jj] = True

        ## evaluate the model
        rownonzero, colnonzero = np.nonzero(outarr['mask'].data)
        lat = outarr.coords['latitude'].values[rownonzero]
        lon = outarr.coords['longitude'].values[colnonzero]
        depths = depth_in_km*np.ones_like(lon)
        if kwargs:
            values = self.evaluate_at_location(lat,lon,depths,parameter, **kwargs)
        else:
            values = self.evaluate_at_location(lat,lon,depths,parameter)
        totarea=0.
        average=0.
        for ii,row in enumerate(rownonzero):
            col = colnonzero[ii]
            outarr['inside'][row,col] = values[ii]
            totarea = totarea + area[row,col].item()
            average = average + values[ii] * area[row,col].item()
        percentarea = totarea/np.sum(area).item()*100.
        average = average/totarea

        ## evaluate values outside polygon
        if outside:
            rowzero, colzero = np.nonzero(np.logical_not(outarr['mask']).data)
            lat = outarr.coords['latitude'].values[rowzero]
            lon = outarr.coords['longitude'].values[colzero]
            depths = depth_in_km*np.ones_like(lon)
            if kwargs:
                values = self.evaluate_at_location(lat,lon,depths,parameter, **kwargs)
            else:
                values = self.evaluate_at_location(lat,lon,depths,parameter)
            totarea_outside=0.
            average_outside=0.
            for ii,row in enumerate(rowzero):
                col = colzero[ii]
                outarr['outside'][row,col] = values[ii]
                totarea_outside = totarea_outside + area[row,col].item()
                average_outside = average_outside + values[ii] * area[row,col].item()
            average_outside = average_outside/totarea_outside
        if outside:
            return percentarea,outarr,average,average_outside
        else:
            return percentarea,outarr,average

    def check_unit(self,units,parameter,resolution):
        if units == None: return
        resolution = tools.convert2nparray(resolution,int2float=False)
        checks=True

        for _,res in enumerate(resolution):
            if 'attrs' not in self[res].keys():
                checks=False
            else:
                if parameter not in self[res]['attrs'].keys():
                    checks=False
                else:
                    if 'unit' not in self[res]['attrs'][parameter].keys(): checks=False
            if not checks: raise ValueError('decoding units not possible without attrs unit initialized for '+parameter+' from model file and attributes.ini')
        if units not in ['absolute','default']:
            res0_unit = self[0]['attrs'][parameter]['unit']
            try:
                junk = constants.ureg(res0_unit).to(units)
            except DimensionalityError:
                raise ValueError("Units can only be None, absolute, default or something that can be converted from default in model instance")
        for index,res in enumerate(resolution):
            if index == 0:
                unit = self[res]['attrs'][parameter]['unit']
            else:
                if unit != self[res]['attrs'][parameter]['unit']:
                    raise AssertionError('units of parameter '+parameter+' need to be same for all resolutions, not '+unit+' and '+self[res]['attrs'][parameter]['unit'])

            if units == 'absolute':
                required = ['absolute_unit','refvalue','start_depths','end_depths'] if self._type == 'netcdf4' else ['absolute_unit']
                for field in required:
                    if field not in self[res]['attrs'][parameter].keys():
                        raise AssertionError(field+' of parameter '+parameter+' need to be provided for resolution '+str(res)+' if units argument is '+units)

                # special checks for netcdf
                if self._type == 'netcdf4':
                    if index == 0:
                        refvalue = self[res]['attrs'][parameter]['refvalue']
                        start_depths = self[res]['attrs'][parameter]['start_depths']
                        end_depths = self[res]['attrs'][parameter]['end_depths']
                    else:
                        if np.any(refvalue != self[res]['attrs'][parameter]['refvalue']): raise AssertionError('refvalue of parameter ',refvalue,' need to be same for all resolutions, not ',refvalue,' and ',self[res]['attrs'][parameter]['refvalue'])
                        if self._type == 'netcdf4':
                            if np.any(start_depths != self[res]['attrs'][parameter]['start_depths']): raise AssertionError('start_depths of parameter '+parameter+' need to be same for all resolutions, not ',start_depths,' and ',self[res]['attrs'][parameter]['start_depths'])
                            if np.any(end_depths != self[res]['attrs'][parameter]['end_depths']): raise AssertionError('end_depths of parameter '+parameter+' need to be same for all resolutions, not ',end_depths,' and ',self[res]['attrs'][parameter]['end_depths'])

    def evaluate_unit(self,parameter,values,units,depth_in_km=None,add_reference=True,resolution=0):
        # change units based on options
        if units is None: return values
        if depth_in_km is None and units == 'absolute': raise ValueError('absolute units do nto work withour depth_in_km')
        # all resolution should have same units so using 0
        unit = self[resolution]['attrs'][parameter]['unit']
        try:
            values_unit=values.to(unit)
        except:
            values_unit = values*constants.ureg(unit)
            warnings.warn('Assuming that the model coefficients  for resolution '+str(resolution)+' are '+unit)

        if units == 'default':
            values = values_unit
        elif units == 'absolute':
            # select reference values
            reference = self.get_reference(parameter,depth_in_km,resolution)
            reference = reference.reshape(values.shape,order='C')
            #reference = refvalue[index].reshape((ndep,nlat),order='C')
            # 1 as reference value is added
            values_fraction = 1+values_unit.to('fraction') if add_reference else values_unit.to('fraction')
            values = np.multiply(values_fraction.magnitude,reference)
        else:
            values = values_unit.to(units)
        return values

    def get_reference(self,parameter,depth_in_km,resolution=0,interpolation='linear',dbs_path=tools.get_filedir()):

        # convert to numpy arrays
        depth_in_km = tools.convert2nparray(depth_in_km)

        # try reading reference model first
        dbs_path=tools.get_fullpath(dbs_path)
        ref1Dfile = dbs_path+'/'+self[resolution]['refmodel']
        if os.path.isfile(ref1Dfile):
            ref1d = Reference1D(ref1Dfile)
            values1D = ref1d.evaluate_at_depth(depth_in_km,parameter,interpolation=interpolation)
            return values1D

        # Try reading as netcdf
        if self._type == 'netcdf4':
            values1D = np.zeros(depth_in_km.size)
            refvalue = self[resolution]['attrs'][parameter]['refvalue']
            start_depths = self[resolution]['attrs'][parameter]['start_depths']
            end_depths = self[resolution]['attrs'][parameter]['end_depths']
            absolute_unit = self[resolution]['attrs'][parameter]['absolute_unit']
            interpolant = self[resolution]['interpolant']
            if interpolant.lower() == 'nearest':
                index = tools.ifwithindepth(start_depths,end_depths,depth_in_km)
                nonzeros = np.where(index>=0)[0]
                values1D[nonzeros]=refvalue[index[nonzeros]]
            elif interpolant.lower() == 'smooth':
                mid_depths = (start_depths+end_depths)/2.
                # Sort depth to make drspln work
                x = mid_depths;xind = x.argsort();x = x[xind]
                x = np.array(x.tolist(), order = 'F')
                y = refvalue; y = y[xind]
                y = np.array(y.tolist(), order = 'F')
                (q,wk)  =  drspln(1,len(x),x,y)
                for ii,depwant in enumerate(depth_in_km):
                    values1D[ii]  =  drsple(1,len(x),x,y,q,depwant)
            else:
                raise KeyError('interpolant type '+interpolant+' cannot be deciphered for get_reference query')
            return values1D*constants.ureg(absolute_unit)
        else:
            raise IOError('Could not fill some reference values as the 1D reference model file could not be read as Reference1D instance : '+ref1Dfile)

    def evaluate_slice(self,parameter,data=None,grid=10.,depth_in_km=None,**kwargs):
        """
        checks whether the data input is a DataArray and the coordinates are compatible

        Parameters
        ----------

        parameter: 2D or 3D parameter

        data : xarray Dataset or a dictionary with values

        grid : size of pixel in degrees, ignore if data is Dataset

        depth_in_km: list of depths in km for a 3D parameter

        Returns
        -------

        data: new or updated xarray Dataset

        """

        # temporary variables
        data_vars={} if data is None else data

        # loop over variables for xarray Dataset
        if isinstance(data, xr.Dataset):
            latitude = data.coords['latitude'].values
            longitude = data.coords['longitude'].values
            old=True
        elif data is None or isinstance(data, dict):
            old=False
            # prepare grid for evaluation
            if isinstance(grid,(int,np.int64,float)):
                latitude = np.arange(-90+grid/2., 90,grid)
                longitude = np.arange(0+grid/2., 360,grid)
                meshgrid =  True

                # check pixel width
                pixel_array= grid*np.ones((len(latitude),len(longitude)))
                if 'pix_width' in data_vars.keys():
                    if data_vars['pix_width'] != pixel_array: raise ValueError('incompatible data_vars')
                else:
                    data_vars['pix_width']=(('latitude', 'longitude'),pixel_array)
            elif isinstance(grid,string_types):
                profiles._infile = grid
                meshgrid =  False
                raise NotImplementedError('needs to be implemented soon')
        else:
            raise ValueError('input data need to be a data_vars dictionary with a grid size or an existing data array or set')

        # check  that incompatible fields are input
        for key in ['grid','interpolated']:
            if key in kwargs.keys(): raise ValueError(key+' cannot be passed as an argument to evaluate_slice')

        # fill values
        if kwargs:
            value = self.evaluate_at_location(parameter,latitude,longitude,depth_in_km,grid=True,interpolated=False,**kwargs)
        else:
            value = self.evaluate_at_location(parameter,latitude,longitude,depth_in_km,grid=True,interpolated=False)

        # store as xarray
        if depth_in_km is None:
            fields = ('latitude', 'longitude')
            coords = dict(zip(fields,(latitude,longitude)))
            data_vars[parameter]=(fields,value[0].magnitude)
        else:
            fields = ('depth','latitude', 'longitude')
            coords = dict(zip(fields,(depth_in_km,latitude,longitude)))
            data_vars[parameter]=(fields,value.magnitude)

        return data_vars if old else xr.Dataset(data_vars = data_vars,coords = coords)


    def evaluate_at_location(self,parameter,latitude,longitude,depth_in_km=None,resolution=0,realization=0,grid=False,interpolated=False,tree=None,nearest=None,units=None,add_reference=True,dbs_path=tools.get_filedir()):
        """
        Evaluate the mode at a location (latitude, longitude,depth)

        Input Parameters:
        ----------------

        parameter: variable whose value will be returned

        depth_in_km: depth where values are being queried at (km)

        interpolated: If True, use KDTree from a predefined grid. If False, evaluated
                      exactly using kernel_set instance.

        nearest: number of points to interpolate. If None, defaults to values based on
                    self[resolution]['interpolant']

        grid: make a grid by unraveling (depth_in_km,latitude,longitude)

        units: if 'absolute' try converting to absolute units of the parameter.
               If 'default', give the default units from model instance.
               If None, simply return the sparse array

        add_reference: add the reference paramter value to perturbation, only valid
                        if units queried is 'absolute'

        """
        if self._name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")

        # convert to numpy arrays
        latitude = tools.convert2nparray(latitude)
        longitude = tools.convert2nparray(longitude)
        resolution = tools.convert2nparray(resolution,int2float=False)
        tree_provided = False if tree == None else True
        nlat = len(latitude)
        nlon = len(longitude)
        if np.all(depth_in_km==None):
            ndep = 1
        else:
            depth_in_km = tools.convert2nparray(depth_in_km)
            ndep = len(depth_in_km)

        #checks
        self.check_unit(units,parameter,resolution)
        if grid:
            nrows = ndep*nlat*nlon
            longitude,latitude = np.meshgrid(longitude,latitude)
            longitude = longitude.ravel()
            latitude = latitude.ravel()
        else:
            check1 = (len(latitude) == len(longitude) == len(depth_in_km)) if depth_in_km is not None else (len(latitude) == len(longitude))
            if not check1: raise ValueError("latitude, longitude or depth_in_km should be of same length if not making grid = False")
            check2 = (latitude.ndim == longitude.ndim == depth_in_km.ndim == 1) if depth_in_km is not None else  (latitude.ndim == longitude.ndim)
            if not check2: raise ValueError("latitude, longitude or depth_in_km should be one-dimensional arrays")
            nrows = nlat

        #compute values for each resolution
        sumvalue = np.zeros((ndep,nlat,nlon),order='C') if grid else np.zeros(nrows)
        for index,res in enumerate(resolution):
            # get the model array
            if parameter in self[resolution]['varstr']:
                mapping={'variables': [parameter], 'constants': [1.0], 'symbols': []}
            else:
                try:
                    mapping = self[resolution]['attrs'][parameter]['mapping']
                    outvar=[];outsym=[];outcons=[]
                    for ivar,variable in enumerate(mapping['variables']):
                        findvar = tools.convert2nparray(self[res]['kernel_set'].search_radial(variable))
                        for itemp,var in enumerate(findvar):
                            outvar.append(var)
                            outcons.append(mapping['constants'][ivar])
                            if ivar > 0:
                                outsym.append(mapping['symbols'][ivar-1])
                            else:
                                if itemp>0: outsym.append('+')
                    mapping['variables']=outvar;mapping['constants']=outcons
                    mapping['symbols']=outsym
                except:
                    raise KeyError('no parameter mapping or entry found for '+parameter+' in resolution '+str(res)+'. Only following queries allowed : '+str(self[resolution]['varstr']))

            # No wloop over mapping and add/subtract as required
            for ipar,param in enumerate(mapping['variables']):
                modelarr = self.coeff2modelarr(resolution=res,realization=realization, parameter=param)

                if not interpolated:
                    # get the projection matrix
                    project = self.get_projection(latitude=latitude,longitude=longitude,depth_in_km=depth_in_km,parameter=param,resolution=res,grid=grid)
                    predsparse = project['matrix']*modelarr

                    # update locations based on grid
                    depth_tmp= None
                    if grid and ipar == 0:
                        if depth_in_km is not None:
                            depth_tmp = np.zeros(nrows)
                            for indx,depth in enumerate(depth_in_km):
                                depth_tmp[indx*nlat*nlon:(indx+1)*nlat*nlon] = depth * np.ones(nlat*nlon)
                ###### interpolant  ####
                else:
                    # decide on the interpolant type based on metadata
                    if nearest == None:
                        interpolant = self[resolution]['interpolant']
                        if interpolant.lower() == 'nearest':
                            nearest=1
                        elif interpolant.lower() == 'smooth':
                            nearest=6
                        else:
                            raise KeyError('interpolant type '+interpolant+' cannot be deciphered for kdtree query')

                    # update locations based on grid
                    if grid:
                        if depth_in_km is None:
                            depth_tmp= None
                        else:
                            depth_tmp = np.zeros(nrows)
                            for indx,depth in enumerate(depth_in_km):
                                depth_tmp[indx*nlat*nlon:(indx+1)*nlat*nlon] = depth * np.ones(nlat*nlon)
                        lat_tmp = np.tile(latitude,len(depth_in_km))
                        lon_tmp = np.tile(longitude,len(depth_in_km))
                    else:
                        depth_tmp = depth_in_km
                        lat_tmp = latitude
                        lon_tmp = longitude

                    # check if the queries is within the bounds of the model
                    checks = self.ifwithinregion(lat_tmp,lon_tmp,depth_tmp,resolution)
                    if np.count_nonzero(checks) == 0:
                        predsparse = sparse.csr_matrix((nrows, 1))
                    else:
                        if not tree_provided:
                            tree = self.buildtree3D(resolution=resolution,dbs_path=dbs_path)
                        self.checktree3D(tree,parameter=param,resolution=resolution)

                        # Get modelarr
                        modelarr = self.coeff2modelarr(resolution=resolution,realization=realization,parameter=param)

                        xlapix = lat_tmp[checks];xlopix = lon_tmp[checks]
                        depths = depth_tmp[checks]
                        #KDtree evaluate
                        print ('... evaluating KDtree ...')
                        temp,_ = tools.querytree3D(tree=tree,latitude=xlapix,longitude=xlopix,radius_in_km= constants.R.to('km').magnitude - depths,values=modelarr,nearest=nearest)
                        print ('... done evaluating KDtree.')
                        if np.all(checks): # if all points within region
                            predsparse = temp
                        else:
                            data = temp.toarray().ravel()
                            row = np.where(checks)[0]
                            col = np.zeros_like(row,dtype=int)
                            predsparse = sparse.csr_matrix((data, (row, col)), shape=(nrows, 1))

                # unravel to get it in the right order
                values = predsparse.toarray().reshape((ndep,nlat,nlon),order='C') if grid else predsparse.toarray().ravel()

                # account for units in model evaluation
                scale=mapping['constants'][ipar]
                if ipar==0:
                    sumvalue+=scale*self.evaluate_unit(param,values,'default',depth_tmp)
                else:
                    if mapping['symbols'][ipar-1]=='+':
                        sumvalue+=scale*self.evaluate_unit(param,values,'default',depth_tmp)
                    elif mapping['symbols'][ipar-1]=='-':
                        sumvalue-=scale*self.evaluate_unit(param,values,'default',depth_tmp)
                    else:
                        raise ValueError('invalid sign '+mapping['symbols'][ipar-1])

        # add backreference values
        if add_reference:
            # first convert to standard units in first resolution
            sumvalue = self.evaluate_unit(parameter,sumvalue,units,depth_tmp,add_reference=True,resolution=0)

        if interpolated and not tree_provided:
            return sumvalue, tree
        else:
            return sumvalue

    def getpixeldepths(self,resolution,parameter):
        typehpar = self[resolution]['typehpar']
        if not len(typehpar) == 1: raise AssertionError('only one type of horizontal parameterization allowed')
        for types in typehpar:
            if not types == 'PIXELS': raise AssertionError('for interpolation with tree3D')
        kernel_set = self[resolution]['kernel_set']
        depths = kernel_set.pixeldepths(parameter)
        return depths

    def coeff2modelarr(self,resolution=0,realization=0,parameter=None):
        """
        Convert the coeffiecient matrix from the file to a sparse model array. Add
        up all the resolutions if a list is provided.

        realization : index of the set of coefficients in an ensemble. Default is 0
        as there is only one set of coefficients when read from a model file.

        resolution : list of resolutions to include the the modelarray

        parameter: parameters to select. Default is to include all available.

        """
        if self._name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")


        # convert to numpy arrays
        resolution = tools.convert2nparray(resolution,int2float=False)

        # Loop over resolution levels, adding coefficients
        for ir,_ in enumerate(resolution):
            try:
                coefficients =self[resolution[ir],realization]['coef']

                #if modelarr is already made, and not parameter is specified, use stored
                recalculate = False
                if parameter == None:
                    try:
                        modelarr = self[resolution[ir],realization]['modelarr']
                    except:
                        recalculate=True
                        Nhopar=len(self[resolution[ir]]['ihorpar'])
                        radselect = np.arange(Nhopar)
                else:
                    recalculate=True
                    # select row indices if specific parameter is required
                    varindex  = np.where(self[resolution[ir]]['varstr']==parameter)[0]+1
                    if len(varindex) != 1: raise AssertionError(str(len(varindex))+' parameters found for '+parameter+'. Only one is allowed')
                    radselect = np.where(self[resolution[ir]]['ivarkern']==varindex[0])[0]

                if recalculate:
                    # if only a single parameterization, simply unravel
                    if len(np.unique(self[resolution]['ihorpar'][radselect])) == 1:
                        modelarr=coefficients.iloc[radselect].to_numpy().ravel()
                    else:
                        for irow,ihorpar in enumerate(self[resolution[ir]]['ihorpar']):
                            # only include if selected
                            icount=0
                            if irow in radselect:
                                ncoef = self[resolution]['ncoefhor'][ihorpar-1]
                                if icount == 1:
                                    modelarr = coefficients.iloc[irow][:ncoef].to_numpy().ravel()
                                else:
                                    modelarr = np.concatenate((modelarr,coefficients.iloc[irow][:ncoef].to_numpy().ravel()))
                    modelarr = sparse.csr_matrix(modelarr).transpose()
                    # only store if all parameters are requested
                    if parameter == None: self[resolution[ir],realization]['modelarr'] = modelarr
            except AttributeError:
                raise ValueError('resolution '+str(resolution[ir])+' and realization '+str(realization)+' not filled up yet.')
        return modelarr

    def readprojbinary(self,lateral_basis):
        """
        Reads Projection matrix created by plot_3dmod_pm.
        lateral_basis can be M362 or pixel1
        """
        if self._name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")

        #read all the bytes to indata
        outfile = self._name+'.'+lateral_basis+'.proj.bin.h5'
        infile = self._name+'.'+lateral_basis+'.proj.bin'

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
                            refstr = struct.unpack('80s',f.read(80))[0].strip().decode('utf-8'); cc = cc+80
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
            model = self._name

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
            model=h5f.attrs['model']; lateral_basis=h5f.attrs['param']

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

    def get_projection(self,parameter,latitude,longitude,depth_in_km = None, resolution = 0,grid=False):
        """
        Get the projection matrix from a lateral basis to another and for particular depths

        depth_in_km : depth in km where the projection matrix is needed.
                      If None, returns the projection matrix for the lat/lon
                      and radial basis as a dirac delta.

        grid: make a grid by repeating (latitude,longitude) by number of depth_in_km

        """
        if self._name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        if latitude is None and longitude is None and depth_in_km is None: raise ValueError("Need to provide atleast one of latitude, longitude or depth_in_km")

        # convert to numpy arrays
        lat = tools.convert2nparray(latitude)
        lon = tools.convert2nparray(longitude)
        if depth_in_km is None:
            dep = None
            nrows = len(lat)
        else:
            dep = tools.convert2nparray(depth_in_km)
            if grid:
                nrows = len(lat)*len(dep)
            else:
                if not (len(lat) == len(lon) == len(dep)): raise ValueError("latitude, longitude or depth_in_km should be of same length if not making grid = False")
                if not (lat.ndim == lon.ndim == dep.ndim == 1): raise ValueError("latitude, longitude or depth_in_km should be one-dimensional arrays")
                nrows = len(lat)

        if not len(latitude)==len(longitude): raise AssertionError('latitude and longitude should be of same length')

        # Get the radial projection file
        kernel = self[resolution]['kernel_set']

        # Write to a dictionary
        # write to a projection file if it does not exist or existing one has different attributes than projection[parameter]
        projection = {}
        projection['kernel'] = kernel.name
        projection['deptharr']=dep
        projection['xlat']=lat
        projection['xlon']=lon
        projection['refstrarr']=parameter

        # loop through parameter and append the projection for each location
        horcof, vercof = kernel.evaluate_bases(parameter,lat,lon,dep)

        # find radial indices of a given physical parameter
        findrad = kernel.find_radial(parameter)

        # convert vercof to sparse for easy multiplication
        vercof = sparse.csr_matrix(vercof) if vercof is not None else np.ones((1,len(findrad)))

        # initialize the projection matrix
        indend_all = kernel.metadata['ncoefcum'][findrad['index'][-1]]-1
        indstart_all = 0 if findrad['index'][0] == 0 else kernel.metadata['ncoefcum'][findrad['index'][0]-1]
        ncol = indend_all - indstart_all + 1
        proj = sparse.lil_matrix((nrows,ncol), dtype=np.float64)

        # difference projection matrix based on whether grid is asked or not
        # loop over all depths
        for idep in range(1 if dep is None else len(dep)):
            # loop over all radial kernels that belong to this parameter and add up
            for ii in np.arange(len(findrad)):
                if vercof[idep,ii] != 0.:
                    row_start = idep*len(lat) if grid else idep
                    row_end = (idep+1)*len(lat) if grid else idep+1
                    column_start = ii*horcof.shape[1]
                    column_end = (ii+1)*horcof.shape[1]
                    proj[row_start:row_end,column_start:column_end] = horcof*vercof[idep,ii] if grid else (horcof*vercof[idep,ii])[idep]

        # store the projection matrix
        projection['matrix'] = proj.tocsr()

        # store the lateral and radial parameter of this variable
        radial_type = []
        for rker in kernel.data['radial_basis'][parameter]: radial_type.append(rker.type)
        if len(np.unique(radial_type)) != 1: raise AssertionError('more than one radial parameterization for  '+parameter)
        projection['radial_basis']=radial_type[0]
        projection['radial_metadata']=kernel.data['radial_basis'][parameter][0].metadata

        # same check for lateral parameterization
        selfmeta = self[resolution]
        selfkernel = selfmeta['kernel_set']
        ivarfind =np.where(selfmeta['varstr']==parameter)[0]
        if not len(ivarfind) == 1: raise AssertionError('only one parameter can be selected in eval_kernel_set')
        findvar = selfmeta['varstr'][ivarfind[0]]
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        findrad = np.array([(ii, selfmeta['desckern'][ii]) for ii in np.arange(len(selfmeta['ivarkern'])) if ivarfind[0]+1 == selfmeta['ivarkern'][ii]],dtype=dt)

        if len(np.unique(selfmeta['ihorpar'][findrad['index']])) != 1: raise AssertionError('more than one lateral parameterization for  '+parameter)
        ihorpar = selfmeta['ihorpar'][findrad['index']][0]
        projection['lateral_metadata']= selfkernel.data['lateral_basis'][ihorpar-1].metadata
        projection['lateral_basis']= selfkernel.data['lateral_basis'][ihorpar-1].type

        return projection

    def projslices(self,projection,variable,depth,resolution=0,realization=0):
        """
        Projection matrix multiplied by model ensemble. Choses the nearest depth available for the projection.
        """
        if self._name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")

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

    def checktree3D(self,tree,parameter,resolution=0):

        # check if it is pixel based
        typehpar = self[resolution]['typehpar']
        if not len(typehpar) == 1: raise AssertionError('only one type of horizontal parameterization allowed')
        for types in typehpar:
            if not types == 'PIXELS': raise AssertionError('for interpolation with tree3D')

        # find the radial kernels
        varindex  = np.where(self[resolution]['varstr']==parameter)[0]+1
        if len(varindex) != 1: raise AssertionError(str(len(varindex))+' parameters found for '+parameter+'. Only one is allowed')
        radselect = np.where(self[resolution]['ivarkern']==varindex[0])[0]
        # if only a single parameterization, simply unravel
        if len(np.unique(self[resolution]['ihorpar'][radselect])) != 1:
            raise ValueError('cannot deal with more than one horizontal parameterization')

        # get the pixel configuration
        depths = self.getpixeldepths(resolution,parameter)
        if np.any(depths == None): raise ValueError('depths not found for interpolation for variable '+parameter+' in target file.')
        #get full path
        xlapix = self[resolution]['xlapix'][0]
        xlopix = self[resolution]['xlopix'][0]
        nlat = len(xlapix)
        ndep = len(depths)
        depth_in_km = np.zeros(nlat*ndep)
        for indx,depth in enumerate(depths):
            depth_in_km[indx*nlat:(indx+1)*nlat] = depth * np.ones(nlat)
        xlapix = np.tile(xlapix,len(depths))
        xlopix = np.tile(xlopix,len(depths))
        radius_in_km = constants.R.to('km').magnitude - depth_in_km

        if tree.n != len(xlopix):
            raise AssertionError('tree (knots = '+str(tree.n)+') not compatible with the self instance for variable '+parameter+' (knots = '+str(len(xlopix))+'). Perhaps fewer than necessary depths have been provided; '+str(ndep)+' depths are available now.')

        # tree stucture is store in spherical coordinates
        rlatlon = np.column_stack((radius_in_km.flatten(),xlapix.flatten(), xlopix.flatten()))
        xyz = spher2cart(rlatlon)
        if not np.allclose(tree.data,xyz): raise ValueError('tree is not compatible with the matrix of the variable '+parameter)

        return xlapix,xlopix,depth_in_km


    def reparameterize(self, model3d,resolution=0,realization=0,interpolated=False,tree=None,nearest=1, dbs_path=tools.get_filedir()):
        """
        Inverts for new coefficients in self.data from the coefficients in model3d class

        write : output the reparamterized model as a text file

        interpolated: assume that it won't be a simple interpolation
        """
        if type(model3d).__name__ != 'Model3D': raise ValueError('input model class should be a Model3D instance')

        # Get the radial projection file
        selfmeta = self[resolution]
        newmeta = model3d[resolution]
        selfkernel = selfmeta['kernel_set']
        newkernel = newmeta['kernel_set']
        tree_provided = False if tree == None else True

        ####################### If both models are pixel based then interpolate###########
        check1 = selfkernel.metadata['typehpar'] == newkernel.metadata['typehpar'] == 'PIXELS'
        # they should also have the parameter in common
        check2 = sorted(np.intersect1d(selfkernel.metadata['varstr'],newkernel.metadata['varstr'])) == sorted(selfkernel.metadata['varstr'])

        if len(check1) == 1 and check1[0] and check2 and interpolated:
            # storage of kernel descriptions
            desckern =  np.array([],dtype=newmeta['desckern'].dtype)
            ivarkern=np.array([],dtype=int); ihorpar=np.array([],dtype=int);
            ncoeff_layer=np.array([],dtype=int)
            for indx, variable in enumerate(selfmeta['varstr']):
                if not tree_provided:
                    tree = self.buildtree3D(resolution=resolution,dbs_path=dbs_path)
                # check if tree is compatible with modelarr
                self.checktree3D(tree,parameter=variable,resolution=resolution)

                # Get modelarr
                modelarr = self.coeff2modelarr(resolution=resolution,realization=realization,parameter=variable)

                # queried locations
                depth_shared = model3d.getpixeldepths(resolution,variable)
                if np.any(depth_shared == None): raise ValueError('depths not found for interpolation for variable '+variable+' in target file.')
                xlapix = newmeta['xlapix'][0]
                xlopix = newmeta['xlopix'][0]
                nlat = len(xlapix)
                nlon = len(xlapix)
                ndep = len(depth_shared)
                depth_in_km = np.zeros(nlat*ndep)
                for indx1,depth in enumerate(depth_shared):
                    depth_in_km[indx1*nlat:(indx1+1)*nlat] = depth * np.ones(nlat)
                xlapix = np.tile(xlapix,len(depth_shared))
                xlopix = np.tile(xlopix,len(depth_shared))

                # check if the queries is within the bounds of the model
                checks = self.ifwithinregion(xlapix,xlopix,depth_in_km,resolution)
                if np.count_nonzero(checks) == 0:
                    final = sparse.csr_matrix((ndep*nlat, 1))
                else:
                    xlapix = xlapix[checks];xlopix = xlopix[checks]
                    depth_in_km = depth_in_km[checks]

                    #KDtree evaluate
                    print ('... evaluating KDtree ...')
                    temp,_ = tools.querytree3D(tree=tree,latitude=xlapix,longitude=xlopix,radius_in_km= constants.R.to('km').magnitude - depth_in_km,values=modelarr,nearest=nearest)
                    print ('... done evaluating KDtree.')
                    if np.all(checks): # if all points within region
                        final = temp
                    else:
                        data = temp.toarray().ravel()
                        row = np.where(checks)[0]
                        col = np.zeros_like(row,dtype=int)
                        final = sparse.csr_matrix((data, (row, col)), shape=(ndep*nlat, 1))

                # now unwrap to coef
                shape = [ndep,nlat]
                if indx == 0:
                    values = final.reshape(shape,order='C')
                else:
                    values = sparse.vstack([values,final.reshape(shape,order='C')])

                # update the kernel descriptions
                varind = np.where(newmeta['varstr']==variable)[0][0]+1
                kerind = np.where(newmeta['ivarkern']==varind)

                ivarkern = np.concatenate((ivarkern, np.ones(len(newmeta['ivarkern'][kerind]),dtype=int)*(indx+1)))
                ihorpar = np.concatenate((ihorpar,newmeta['ihorpar'][kerind]))
                desckern = np.concatenate((desckern,newmeta['desckern'][kerind]))
                ncoeff_layer = np.concatenate((ncoeff_layer, newmeta['ncoefhor'][newmeta['ihorpar'][kerind]-1]))

                #update the variable attributes, delete some values as they will be recalculated
                for field in newmeta['attrs'][variable].keys():
                    if not field in ['refvalue','average']:
                        selfmeta['attrs'][variable][field] = newmeta['attrs'][variable][field]
                    else:
                        if field in selfmeta['attrs'][variable]:
                            del selfmeta['attrs'][variable][field]

            # copy over the correct metadata
            selfmeta['nmodkern'] = len(selfmeta['desckern'])
            selfmeta['ncoefcum'] = np.cumsum(ncoeff_layer)
            selfmeta['ivarkern'] = ivarkern
            selfmeta['ihorpar'] = ihorpar
            selfmeta['desckern'] = desckern
            for field in ['hsplfile','xsipix', 'xlapix','ncoefhor', 'xlopix','kerstr','geospatial_lat_resolution','geospatial_lon_resolution']:
                selfmeta[field] = newmeta[field]
            try:
                selfmeta['kernel_set'] = Kernel_set(selfmeta.copy())
            except:
                warnings.warn('Warning: kernel_set could not initialized for '+str(resolution))
                pass

            # make a model3D instance and store coef panda dataframe
            reparam = copy.deepcopy(self)
            reparam[resolution] = selfmeta
            reparam[resolution,realization]['coef'] = pd.DataFrame(values.toarray())
            reparam._infile = selfmeta['name']+'.'+selfmeta['kerstr']+'.avni.nc4'

        ####################### Invert the coefficients    ##############################
        else:
            if interpolated: warnings.warn('Could not interpolate. Inverting the coefficients...')
            # get the projection matrix for each variable in self
            dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
            for variable in selfmeta['varstr']:
                ivarfind =np.where(selfmeta['varstr']==variable)[0]
                if not len(ivarfind) == 1: raise AssertionError('only one parameter can be selected in eval_kernel_set')
                findvar = selfmeta['varstr'][ivarfind[0]]
                findrad = np.array([(ii, selfmeta['desckern'][ii]) for ii in np.arange(len(selfmeta['ivarkern'])) if ivarfind[0]+1 == selfmeta['ivarkern'][ii]],dtype=dt)

                # Check if the selected radial kernels are boxcar in self
                for rker in selfkernel.data['radial_basis'][variable]:
                    if rker.type != 'boxcar': raise AssertionError('radial kernel is not boxcar for '+variable)
                # same check for lateral parameterization
                for hpar in selfmeta['ihorpar'][findrad['index']]:
                    if selfkernel.data['lateral_basis'][hpar-1].type != 'PIXELS': raise AssertionError('lateral kernel is not PIXELS for '+variable)

                # find the corresponding radial kernels in newmeta
                ivarfind2 =np.where(newmeta['varstr']==variable)[0]
                if len(ivarfind2) != 1:
                    ifselected = [False for x in range(len(newmeta['varstr']))]
                    ivarfind2 = []
                    ifdone = False
                    while not ifdone:
                        stringout = ''
                        for ii in range(len(newmeta['varstr'])): stringout=stringout+str(ii)+'. '+newmeta['varstr'][ii]+' '+str(ifselected[ii] if ifselected[ii] else '')+'\n'
                        print('')
                        print(stringout)
                        try:
                            x = int(input('Warning: no unique corresponding variable found for '+variable+'. Select one index from above to assign parameterization:'))
                            ifselected[x] = True
                            ivarfind2.append(x)
                        except (ValueError,EOFError):
                            if len(ivarfind2) != 0: ifdone = True

                # loop over selected variables
                for ivar in ivarfind2:
                    findvar2 = newmeta['varstr'][ivar]
                    findrad2 = np.array([(ii, newmeta['desckern'][ii]) for ii in np.arange(len(newmeta['ivarkern'])) if newmeta['ivarkern'][ii] == ivar+1],dtype=dt)

                    # calculate the projection matrix for all locations
                    longitude = selfmeta['xlopix'][0]; latitude = selfmeta['xlapix'][0]
                    radialinfo = selfkernel.data['radial_basis'][variable][0].metadata
                    depth_in_km = np.average(np.array([radialinfo['depthtop'],radialinfo['depthbottom']]),axis=0)
                    # loop over depths and append the projection matrices
                    proj = model3d.get_projection(findvar2,latitude,longitude,depth_in_km,resolution=resolution)

                    # get areas for the grid and multiply the proj

                    # invert for the coefficients, if not invertible perform L curve inversion
                    # find optimal fraction of max(GTG_diag) to use for damping by
                    # inflection point

                    # select the sub-array with non-zero values to invert
                    start = newmeta['ncoefcum'][findrad2['index'][0]-1] if findrad2['index'][0] > 0 else 0
                    end = newmeta['ncoefcum'][findrad2['index'][-1]]

                    GTG= proj['matrix'].T*proj['matrix']
                    GTG_inv = linalg.inv(GTG[start:end,start:end].todense())
                    #d = GTG_inv * values
                    pdb.set_trace()

        if interpolated and not tree_provided:
            return reparam, tree
        else:
            return reparam

    def to_profiles(self,grid=10.,type='pixel',resolution=0,realization=0,model_dir='.',interpolated=True):
        """
        converts a model3d class to profiles class

        grid: either a grid size of a grid file

        interpolant: interpolant type either in pixel or nearest

        model_dir: directory to find reference 1D model

        """
        if type not in ['pixel','nearest']:
            raise ValueError('only pixel or nearest options allowed for type')

        # check if we can read 1D model
        ref1Dfile = self[resolution]['refmodel']
        if os.path.isfile(ref1Dfile):
            try: # try reading the 1D file in card format
                ref1d = Reference1D(ref1Dfile)
                depth_in_km = (ref1d.data['depth'].pint.to('km')).data
                ndep = len(depth_in_km)
            except:
                raise IOError('Could not fill some reference values as the 1D reference model file could not be read as Reference1D instance : '+ref1Dfile)
        else:
            raise IOError('Could not fill some reference values as the 1D reference model file could not be found : '+ref1Dfile)

        tree = None
        # Operations between PintArrays of different unit registry will not work.
        # We can change the unit registry that will be used in creating new
        # PintArrays to prevent this issue.
        pint.PintType.ureg = constants.ureg
        PA_ = pint.PintArray

        # loop over reference1d fields and update the corresponding fields
        common_fields = np.intersect1d(self[resolution]['parameters'],ref1d.metadata['parameters'])
        missing_fields = sorted(set(ref1d.metadata['parameters'])-set(common_fields))

        # copy the 1D model and drop missing fields not specified by 3d model
        local=copy.deepcopy(ref1d)
        local.data=local.data.drop(missing_fields,axis=1)
        local.metadata['parameters'] = common_fields
        index = []
        for field in common_fields:index.append(ref1d.metadata['parameters'].index(field))
        local.metadata['units'] = (np.array(ref1d.metadata['units'])[index]).tolist()
        for field in ['delta','average','contrast']:
            local.metadata['discontinuities'][field] = local.metadata['discontinuities'][field].drop(missing_fields,axis=1)

        # loop over variables for xarray Dataset
        data_vars={};local_profiles = {}
        profiles = Profiles()
        profiles._name = self._name
        profiles._interpolant = type

        if isinstance(grid,(int,np.int64,float)):
            profiles._infile = str(grid)+'X'+str(grid)
            latitude = np.arange(-90+grid/2., 90,grid)
            longitude = np.arange(0+grid/2., 360,grid)
            meshgrid =  True
            nlat = len(latitude)
            nlon = len(longitude)
            nprofiles = nlat*nlon
            # index column by longitude
            profindx = np.arange(nprofiles).reshape((nlat,nlon),order='C')
            data_vars['index']=(('latitude', 'longitude'), profindx)
            # deep copy will not modify values of original profile
            # get the grid sizes stored
            pixel_array= grid*np.ones((len(latitude),len(longitude)))
            data_vars['pix_width']=(('latitude', 'longitude'), pixel_array)
        elif isinstance(grid,string_types):
            profiles._infile = grid
            meshgrid =  False
            raise NotImplementedError('needs to be implemented soon')

        for ii in np.arange(nprofiles): local_profiles[ii]=copy.deepcopy(local)
        # get the grid sizes stored
        profiles.data['grid'] = xr.Dataset(data_vars = data_vars,
                coords = {'latitude':latitude,'longitude':longitude})

        tree = None # create the tree the first time arounf
        for indx,parameter in enumerate(common_fields):
            print('... evaluating 3D perturbations for parameter # '+str(indx+1)+' / '+str(len(common_fields))+' '+parameter)
            # Get 3D perturbations, do not add approximate reference stored ins file
            if interpolated:
                if tree == None:
                    values3D, tree = self.evaluate_at_location(longitude=longitude,latitude=latitude,depth_in_km=depth_in_km,parameter=parameter,resolution=resolution,realization=realization,grid=meshgrid,interpolated=True,units='absolute',add_reference=False)
                else:
                    values3D = self.evaluate_at_location(longitude=longitude,latitude=latitude,depth_in_km=depth_in_km,parameter=parameter,resolution=resolution,realization=realization,grid=meshgrid,interpolated=True,tree=tree,units='absolute',add_reference=False)

            # calculate explicitly
            else:
                values3D = self.evaluate_at_location(longitude=longitude,latitude=latitude,depth_in_km=depth_in_km,parameter=parameter,resolution=resolution,realization=realization,grid=meshgrid,interpolated=False,units='absolute',add_reference=False)

            # Get 1D values
            values1D = ref1d.evaluate_at_depth(depth_in_km,parameter)

            # Add the 1D values
            for idep in np.arange(ndep): values3D[idep,:,:] = values3D[idep,:] + values1D[idep]

            # loop over profiles and update values
            if meshgrid:
                done = []
                for ilat in np.arange(nlat):
                    for ilon in np.arange(nlon):
                        index = profiles.data['grid']['index'][ilat,ilon].item()
                        #if not evaluated before
                        if index not in done:
                            target_unit = local_profiles[index].data[parameter].pint.units
                            local_profiles[index].data[parameter][:] = values3D[:,ilat,ilon].to(target_unit)
                            done.append(index)
                            # update name
                            if indx == 0: local_profiles[index]._name = self._name+'_'+ \
                                            local_profiles[index]._name + \
                                            '_profile#' + str(index)
            else:
                raise NotImplementedError('needs to be implemented soon')

        ## Update these values
        for key in self[resolution].keys():
            value = self[resolution][key]
            if isinstance(value, (list,tuple,np.ndarray,bool,float,int,np.int64,string_types)):
                profiles.metadata[key] = value
        profiles.data['profiles'] = local_profiles
        return profiles

    def get_resolution(self,rescovfile=None,LU2symmetric=True,resolution=0, realization=0):
        """
        Reads Resolution or Covariance matrix created by invwdata_pm64 with option -r.
        R=inv(ATA+DTD)ATA and the name of file is typically outmodel.Resolution.bin
        First read the matrix and model file, perform checks and create a sparse matrix.
        """
        # Read metadata from the 3D model
        refmodel1 = self[resolution]['refmodel'];kerstr1 = self[resolution]['kerstr']
        modelarr = self.coeff2modelarr(resolution,realization)

        # default location of resolution file
        if rescovfile is None: rescovfile = self._infile+'.Resolution.bin'
        outfile=rescovfile+'.npz'

        if (not os.path.isfile(outfile)):
            if (not os.path.isfile(rescovfile)): raise IOError("Filename ("+rescovfile+") does not exist")
            # Read the binary covariance matrix, check for metadata
            refmodel2, kerstr2, ntot, indexrad1, indexrad2, indexhor1, indexhor2, out = \
                readResCov(rescovfile,onlymetadata=True)

            # Checks
            if refmodel1 != refmodel2: warnings.warn(" refmodel in model3d instance : "+refmodel1+" not equal to the one in the binary file: "+refmodel2)
            if kerstr1 != kerstr2: raise ValueError("kerstr: "+kerstr1+" not equal to the one in the binary file: "+kerstr2)
            #if ntot != modelarr.size: raise ValueError("Number of variables: ",str(ntot)," not equal to ",str(modelarr.size))

            # Read the binary covariance matrix, check for metadata
            refmodel2, kerstr2, ntot, indexrad1, indexrad2, indexhor1, indexhor2, out = \
                readResCov(rescovfile,onlymetadata=False)

            # Create the sparse matrix
            ncoefcum = self[resolution]['ncoefcum'] # cumulative number of coeffients for each radial kernel
            row=np.zeros(len(indexrad2),dtype=int);col=np.zeros(len(indexrad2),dtype=int)
            for ii in np.arange(len(indexrad2)):
                row[ii] = ncoefcum[indexrad2[ii]-2]+(indexhor2[ii]-1) if indexrad2[ii] > 1 else indexhor2[ii]-1
                col[ii] = ncoefcum[indexrad1[ii]-2]+(indexhor1[ii]-1) if indexrad1[ii] > 1 else indexhor1[ii]-1
            outsparse = sparse.csr_matrix((out, (row, col)), shape=(modelarr.shape[0], modelarr.shape[0]))
            sparse.save_npz(outfile,outsparse)
            print(".... written "+outfile)
        else:
            print(".... reading "+outfile)
            outsparse=sparse.load_npz(outfile)

        # Get the full symmetric matrix
        if LU2symmetric: outsparse=getLU2symmetric(outsparse)
        return outsparse,modelarr

    def printsplinefiles(self):
        """
        Prints out the splines knots into a file.

        Parameters
        -----------

        model3d : The model object from read3dmodelfile

        """
        if self._name == None: raise ValueError("No three-dimensional model has been read into this model3d instance yet")
        for ii in np.arange(len(self.metadata)):

            ifindspline=np.where(self[ii]['typehpar']=='SPHERICAL SPLINES')[0]
            for ifind in ifindspline:
                filename=self[ii]['hsplfile'][ifind]
                arr = np.vstack([self[ii]['xlospl'][ifind],self[ii]['xlaspl'][ifind],self[ii]['xraspl'][ifind]]).transpose()
                np.savetxt(filename,arr, fmt='%7.3f %7.3f %7.3f')
                print(".... written "+filename)
        return
#######################################################################################

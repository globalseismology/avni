#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard AVNI format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys,os
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int,list,tuple

import traceback
from configobj import ConfigObj
import warnings
import numpy as np #for numerical analysis
import pdb
import h5py
import xarray as xr
import pandas as pd

####################### IMPORT AVNI LIBRARIES  #######################################
from .reference1d import Reference1D
from .common import readepixfile
from .. import tools
#######################################################################################

# Profiles of 1D Earth models
class Profiles(object):
    '''
    A class for the profiles from a laterally heteroegeneous model
    '''
    #########################       magic       ##########################
    def __init__(self,file=None,**kwargs):
        self.metadata ={}
        self.data = {}
        self._name = None
        self._interpolant = None
        self._infile = None
        if file is not None: self.read(file,**kwargs)

    def __str__(self):
        if self._name is not None:
            output = "%s is a set of profiles read from %s" % (self._name,self._infile)
        else:
            output = "No profiles have been read into this profiles instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._name})'.format(self=self)

    def __add__(self, other):
        raise NotImplementedError('method to add profiles on top of each other. Should use the add method in reference1D')

    def __getitem__(self,key):
        """returns data from a profile with key"""
        return self.data['profiles'][key]

    #########################       decorators       ##########################

    @property
    def name(self):
        return self._name

    @property
    def interpolant(self):
        return self.metadata['interpolant']

    #########################       methods       #############################

    def read(self,file,**kwargs):
        # try setup.cfg on this folder
        self.read_setup(file)

    def read_setup(self,file='setup.cfg'):
        '''
        Try reading a folder containing ascii files for every location on the surface

        Input parameters:
        ----------------

        setup_file: setup file containing metadata for the model

        '''
        if not os.path.isfile(file):
            raise IOError('No configuration file found.'\
                 'Model directory must contain '+file)
        else:
            parser = ConfigObj(file)
        for key in parser['metadata'].keys(): self.metadata[key] = parser['metadata'][key]

        # loop over all files as avni.models.reference1d
        #   save as a list of classes
        epixarr,metadata,comments = readepixfile(self.metadata['index'])
        profiles = {}
        location_indices = np.unique(epixarr['val'])
        for loc in location_indices:
            file_name = self.metadata['folder'] + '/' + self.metadata['prefix'] + str(loc)
            if os.path.isfile(file_name): profiles[loc] = Reference1D(file_name)
        if len(profiles) != len(location_indices): warnings.warn('Only '+str(len(profiles))+' profiles out of '+str(len(location_indices))+' distinct profiles have been found and read')
        self.data['profiles'] = profiles
        self.data['grid'] = epixarr
        self._name = self.metadata['name']
        self._interpolant = self.metadata['interpolant']
        self._infile = file

    def writehdf5(self,outfile = None, overwrite = False):
        """writes profile class to an hdf5 container"""
        if outfile == None: outfile = self._name+'.profiles.h5'
        #raise NotImplementedError('method to add profiles')
        if overwrite:
            hf = h5py.File(outfile, 'w')
        else:
            hf = h5py.File(outfile, 'a')
        g1 = hf.require_group(self._interpolant)
        if self._name != None: g1.attrs['name']=self._name
        if self._infile != None: g1.attrs['infile']=self._infile
        #pdb.set_trace()

        #save metadata
        hf.create_group('metadata')
        for key, item in self.metadata.items():
            if type(item) == list:
                asciiList = [n.encode("ascii","ignore") for n in item]
                hf['metadata'][key] = asciiList
            else:
                try:
                    hf['metadata'][key] = item
                except TypeError:
                    asciiList = [n.encode("ascii","ignore") for n in item]
                    hf['metadata'][key] = asciiList

        #create group for profile data
        data_group = hf.create_group('data')

        #save grid data
        grid_group = data_group.create_group('grid')
        grid_group['latitude'] = self.data['grid'].latitude.data
        grid_group['longitude'] = self.data['grid'].longitude.data
        grid_group['pix_width'] = self.data['grid'].pix_width.data
        grid_group['index'] = self.data['grid'].index.data

        #save profile data
        all_profiles = data_group.create_group('profiles')

        #self.data['profiles'] is a dictionary
        for i in self.data['profiles']:
            pf = self.data['profiles'][i]
            h5_profile = all_profiles.create_group('{}'.format(i)) #create profile group for current index

            #write profile data
            h5_profile.create_group('data')
            for item in pf.data.keys():
                h5_profile['data'][item] = pf.data[item][:].to_numpy(dtype='float32')

            #create profile metadata
            h5_profile.create_group('metadata')

            #discontinuities is a nested dictionary
            disc = h5_profile['metadata'].create_group('discontinuities')
            disc_delta = disc.create_group('delta')
            disc_average = disc.create_group('average')
            disc_contrast = disc.create_group('contrast')

            for key,item in pf.metadata.items():
                print(key)
                if type(item) == list:
                    asciiList = [n.encode("ascii","ignore") for n in item]
                    h5_profile['metadata'][key] = asciiList
                #elif item == None:

                elif key ==  'discontinuities':
                    for item_ in pf.metadata[key].keys():

                        if item_ == 'delta':
                            print('DELTA')
                            disc_delta['depth'] = pf.metadata[key][item_]['depth'].data
                            disc_delta['radius'] = pf.metadata[key][item_]['radius'].data
                            disc_delta['vsv'] = pf.metadata[key][item_]['vsv'].data
                            disc_delta['vsh'] = pf.metadata[key][item_]['vsh'].data
                            disc_delta['vs'] = pf.metadata[key][item_]['vs'].data
                        elif item_ == 'average':
                            disc_average['depth'] = pf.metadata[key][item_]['depth'].data
                            disc_average['radius'] = pf.metadata[key][item_]['radius'].data
                            disc_average['vsv'] = pf.metadata[key][item_]['vsv'].data
                            disc_average['vsh'] = pf.metadata[key][item_]['vsh'].data
                            disc_average['vs'] = pf.metadata[key][item_]['vs'].data
                        elif item_ == 'contrast':
                            disc_contrast['depth'] = pf.metadata[key][item_]['depth'].data
                            disc_contrast['radius'] = pf.metadata[key][item_]['radius'].data
                            disc_contrast['vsv'] = pf.metadata[key][item_]['vsv'].data
                            disc_contrast['vsh'] = pf.metadata[key][item_]['vsh'].data
                            disc_contrast['vs'] = pf.metadata[key][item_]['vs'].data
                        else:
                            disc[item_] = pf.metadata[key][item_]

                        continue

                elif isinstance(item,type(None)):
                    print('{} is NoneType... not being written'.format(key))
                else:
                    try:
                        h5_profile['metadata'][key] = item
                    except TypeError:
                        asciiList = [n.encode("ascii","ignore") for n in item]
                        h5_profile['metadata'][key] = asciiList

#             for key in self.metadata['resolution_'+str(ires)].keys():
#                 keyval = self.metadata['resolution_'+str(ires)][key]
#                 try:
#                     if keyval.dtype.kind in {'U', 'S'}:
#                         keyval = [n.encode('utf8') for n in keyval]
#                         g2.attrs[key]=keyval
#                 except:
#                     if keyval != None and key != 'kernel_set':
#                         g2.attrs[key]=keyval

    def readhdf5(self,hf,interpolant=None,only_metadata = False):
        """
        Reads a standard 3D model file from a hdf5 file

        Input Parameters:
        ----------------

        interpolant: group to load

        only_metadata: do not return the pandas dataframe if True
        """
        f = h5py.File(hf,'r')
        self._name = f['metadata']['name'].value
        self.data['profiles'] = {}

        #read metadata
        for item in f['metadata'].keys():
            self.metadata[item] = f['metadata'][item].value

        #read grid data into xarray
        grid = xr.Dataset()
        grid['latitude'] = f['data']['grid']['latitude'].value
        grid['longitude'] = f['data']['grid']['longitude'].value
        grid['index'] = (('latitude','longitude'),f['data']['grid']['index'].value)
        grid['pix_width'] = (('latitude','longitude'),f['data']['grid']['pix_width'].value)
        self.data['grid'] = grid

        #read individual profiles
        for item in f['data']['profiles']:

            #create an instance of the Reference1D class
            r1d = Reference1D()
            r1d.metadata['discontinuities'] = {}
            r1d._name = '{}_{}_profile#{}'.format(f['metadata']['name'].value,
                                                  f['metadata']['refmodel'].value,
                                                  item)

            #add metadata
            for key_ in f['data']['profiles'][item]['metadata'].keys():
                if key_ == 'discontinuities':
                    for disc_key in f['data']['profiles'][item]['metadata'][key_].keys():

                        if disc_key == 'average':
                            avg_df = pd.DataFrame()
                            avg_df['depth'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['depth'].value
                            avg_df['radius'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['radius'].value
                            avg_df['vsv'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vsv'].value
                            avg_df['vsh'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vsh'].value
                            avg_df['vs'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vs'].value
                            r1d.metadata['discontinuities']['average'] = avg_df

                        elif disc_key == 'contrast':
                            con_df = pd.DataFrame()
                            con_df['depth'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['depth'].value
                            con_df['radius'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['radius'].value
                            con_df['vsv'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vsv'].value
                            con_df['vsh'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vsh'].value
                            con_df['vs'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vs'].value
                            r1d.metadata['discontinuities']['contrast'] = con_df

                        elif disc_key == 'delta':
                            del_df = pd.DataFrame()
                            del_df['depth'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['depth'].value
                            del_df['radius'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['radius'].value
                            del_df['vsv'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vsv'].value
                            del_df['vsh'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vsh'].value
                            del_df['vs'] = f['data']['profiles'][item]['metadata'][key_][disc_key]['vs'].value
                            r1d.metadata['discontinuities']['delta'] = del_df

                        else:
                            r1d.metadata['discontinuities'][disc_key] = f['data']['profiles'][item]['metadata'][key_][disc_key]

                else:
                    r1d.metadata[key_] = f['data']['profiles'][item]['metadata'][key_].value

            #add data (pandas dataframe)
            df = pd.DataFrame()
            df['radius'] = f['data']['profiles'][item]['data']['radius'].value
            df['depth'] = f['data']['profiles'][item]['data']['depth'].value
            df['vsv'] = f['data']['profiles'][item]['data']['vsv'].value
            df['vsh'] = f['data']['profiles'][item]['data']['vsh'].value
            df['vs'] = f['data']['profiles'][item]['data']['vs'].value
            r1d.data = df

            #add Reference1D model to Profiles data dictionary
            self.data['profiles'][int(item)] = r1d

    def find_index(self,latitude,longitude):
        """finds the nearest point in self.data['index']"""
        if self._interpolant == 'pixel':
            #find pixel index from xarray of the lat lon
            indices = self.data['grid']['index'].sel(latitude=latitude, longitude = longitude ,method='nearest')
        elif self._interpolant == 'nearest':
            # use voronoi nearest point
            raise NotImplementedError('method to find nearest point')
        else:
            raise ValueError('only pixel or nearest options allowed for interpolant')
        return indices

    def evaluate_at_location(self,latitude,longitude,depth_in_km,parameter,**kwargs):
        """evaluate the profiles at a particular point within the domain"""
        indices = self.find_index(latitude,longitude)
        if indices.count().item() != 1: raise ValueError('only single location can be queried')
        ref1d = self[indices.item()]
        # evaluate ref1d at this depth and variable from
        if kwargs:
            values = ref1d.evaluate_at_depth(depth_in_km,parameter,**kwargs)
        else:
            values = ref1d.evaluate_at_depth(depth_in_km,parameter,boundary='+',interpolation='linear')
        return values

    def get_profile(self,index,parameters):
        "Use evaluate_at_location at depths defined in reference1D"
        if not isinstance(index, (int,np.int64)): raise KeyError('only integer indices allowed to be queried')
        parameters = tools.convert2nparray(parameters,allowstrings=True)
        select = ['radius','depth']
        for parameter in parameters:
            if parameter not in self[index].metadata['parameters']: raise KeyError(parameter+' not found in the profiles')
            select.append(parameter)
        return self[index].data[select]
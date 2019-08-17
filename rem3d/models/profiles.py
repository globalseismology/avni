#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard REM3D format."""

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

####################### IMPORT REM3D LIBRARIES  #######################################
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

        # loop over all files as rem3d.models.reference1d
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
        raise NotImplementedError('method to add profiles')
        if overwrite:
            hf = h5py.File(outfile, 'w')
        else:
            hf = h5py.File(outfile, 'a')
        g1 = hf.require_group(self._interpolant)
        if self._name != None: g1.attrs['name']=self._name
        if self._infile != None: g1.attrs['infile']=self._infile
        pdb.set_trace()

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
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
        return '{self.__class__.__name__}({self._infile})'.format(self=self)

    def __add__(self, other):
        raise NotImplementedError('method to add profiles on top of each other. Should use the add method in reference1D')
    #########################       decorators       ##########################

    @property
    def name(self):
        return self._name

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
        # use reference1d class and read_bases_coefficients function here
        # go through all files in the folder
        # temp1d = Reference1D(file)

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
        self.data['index'] = epixarr
        self._name = self.metadata['name']
        self._interpolant = self.metadata['interpolant']
        self._infile = file


    def write_to_hdf(self,outfile = None, overwrite = False):
        """writes profile class to an hdf5 container"""
        if outfile == None: outfile = self._name+'.profiles.h5'
        raise NotImplementedError('method to add profiles')
        if overwrite:
            hf = h5py.File(outfile, 'w')
        else:
            hf = h5py.File(outfile, 'a')
        g1 = hf.require_group(self._interpolant)


    def find_index(self,latitude,longitude):
        """finds the nearest point in self.data['index']"""
        indx = None
        if self._interpolant == 'pixel':
            #find pixel index from xarray of the lat lon
            raise NotImplementedError('method to find nearest point')
        elif self._interpolant == 'nearest':
            # use voronoi nearest point
            raise NotImplementedError('method to find nearest point')
        else:
            raise ValueError('only pixel or nearest options allowed for interpolant')
        return indx

    def evaluate_at_location(self,parameter,latitude,longitude,depth_in_km,**kwargs):
        """evaluate the profiles at a particular point within the domain"""
        index = self._find_index(latitude,longitude)
        ref1d = self.data['profiles'][index]
        # evaluate ref1d at this depth and variable from
        if kwargs:
            values = ref1d.evaluate_at_depth(depth_in_km,parameter,**kwargs)
        else:
            values = ref1d.evaluate_at_depth(depth_in_km,parameter,boundary='+',interpolation='linear')
        return values

    def get_profiles(self,parameter,latitude,longitude):
        "Use evaluate_at_location at depths defined in reference1D"

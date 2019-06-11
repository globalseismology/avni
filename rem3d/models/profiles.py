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
    def __init__(self,folder=None):
        self.metadata = None
        self.data = None
        self._name = None
        self._infile = None    
        if file is not None:  success = self._read(file)

    def __str__(self):
        if self._name is not None:
            output = "%s is a set of profiles read from %s" % (self._name,self._infile)
        else:
            output = "No profiles have been read into this profiles instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._infile})'.format(self=self)

    #########################       decorators       ##########################

    @property
    def name(self):
        return self._name
        
    #########################       methods       #############################
        
    def read(self,setup_file='setup.cfg',model_dir='.'):
        '''
        Try reading a folder containing ascii files for every location on the surface 

        Input parameters:
        ----------------

        setup_file: setup file containing metadata for the model

        '''
        # use reference1d class and read_bases_coefficients function here
        # go through all files in the folder
        # temp1d = Reference1D(file)

        cfg_file = model_dir+'/'+setup_file

        if not os.path.isfile(cfg_file):
            raise IOError('No configuration file found.'\
	                 'Model directory must contain '+setup_file)
        else:
            parser = ConfigObj(cfg_file)
        self._folder = parser['metadata']['folder']
        self._index = parser['metadata']['index']
        # loop over all files as rem3d.models.reference1d
        #   save as a list of classes
        epixarr,metadata,comments = readepixfile(self._index)
        profiles = {}
        for loc in np.unique(epixarr['val']):
            file_name = self._folder + '/' + parser['metadata']['prefix'] + str(loc)
            profile = reference1D(file_name)
            profiles[(lat,lon)] = profile
             
    
    def write_to_hdf(self):
        """writes profile class to an hdf5 container"""
        
    def evaluate_at_location(self,latitude, longitude,depth):
        """evaluate the profiles at a particular point within the domain"""
        
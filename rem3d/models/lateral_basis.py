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

####################### IMPORT REM3D LIBRARIES  #######################################
from .. import tools
#######################################################################################
# Horizontal basis parameter class that defines an unique combination of parameters, their radial parameterization and any scaling
# 3D model class
class lateral_basis(object):
    '''
    A class for radial bases that defines a unique combination of parameters,
    their radial parameterization and any scaling that is used.
    '''
    def __init__(self, name, type, metadata = None):
        """
        types : 'epix','ylm','sh','wavelet','slepians'
        """
        self.data = {}
        self.name = name
        self.type = type
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def addtypes(self, names, types):
        """
        types = ['ylm','sh','wavelet','slepians']
        """
        # check the type
        if len(names) !=len(types): raise ValueError("len(names) !=len(types)")
        if isinstance(types,string_types): types = np.array(types)
        for ii in np.arange(types.size):
            # if does not exist already`or not the same as self type
            if types[ii] not in self.proj.keys() and types[ii] != self.metadata['type']:
                self.proj[types[ii]] = {}
                self.proj[types[ii]][to_name[ii]] = {'data':None,'attributes':{}}

    def readprojfile(self,infile):
        """
        Reads a projection matrix file that evaluates the radial bases at various depths.
        """

        if (not os.path.isfile(infile)): raise IOError("Filename ("+infile+") does not exist. Use shell script print_projmatrix to create it.")
        #nbytes = os.path.getsize(infile)

        cc = 0 #initialize byte counter
        ifswp = '' # Assuem that byte order is not be swapped unless elat is absurdly high

        with open(infile, "rb") as f:
            # preliminary metadata
            indata = f.read(4); cc = cc+4 # try to read iflag
            iflag = struct.unpack(ifswp+'i',indata)[0] # Read flag
            if iflag != 1:
                ifswp = '!' # swap endianness from now on
                iflag = struct.unpack(ifswp+'i',indata)[0]
                if iflag != 1: sys.exit("Error: iflag != 1")
            self.metadata['from_type'] = struct.unpack('20s',f.read(20))[0].strip(); cc = cc+20
            pdb.set_trace()
            self.proj['to_type'] = struct.unpack('20s',f.read(20))[0].strip(); cc = cc+20
            #ndp = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
            #npx = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4
            #nhorcum = struct.unpack(ifswp+'i',f.read(4))[0]; cc = cc+4


    def eval_lateral(self,lat,lon,store=False):
        """
        Evaluate radial basis at a depth interpolated from existing projection matrices.
        """
        if self.type == 'SPHERICAL SPLINES':
            horcof = tools.eval_splcon(lat,lon,self.metadata['xlaspl'],self.metadata['xlospl'],self.metadata['xraspl'])
        elif self.type == 'SPHERICAL HARMONICS':
            horcof = tools.eval_ylm(lat,lon, self.metadata['lmaxhor'])
        elif self.type == 'PIXELS':
            horcof = tools.eval_pixel(lat,lon, self.metadata['xlapix'],self.metadata['xlopix'],self.metadata['xsipix'])
        else:
            raise NotImplementedError(self.type+' has not been implemented in lateral_basis.eval_lateral')

        if store:
            self.data['latitude'] = lat
            self.data['longitude'] = lon
            self.data['horcof'] = horcof
        else:
            return horcof

    def project_lateral(self,lateral_basis):
        """
        Project from current horizontal basis to another orthogonal basis
        and return the coefficients.
        """

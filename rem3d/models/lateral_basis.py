#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

import os
import numpy as np #for numerical analysis
import pdb    #for the debugger pdb.set_trace()

####################### IMPORT REM3D LIBRARIES  #######################################
from .. import tools
#######################################################################################
# Horizontal basis parameter class that defines an unique combination of parameters, their lateral parameterization

class Lateral_basis(object):
    '''
    A class for radial bases that defines a unique combination of parameters,
    their radial parameterization and any scaling that is used.
    '''
    #########################       magic       ##########################

    def __init__(self, name, types, metadata = None):
        """
        types : 'epix','ylm','sh','wavelet','slepians'
        """
        self.data = {}
        self.name = name
        self.type = types
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def __eq__(self, other):

        # convert to array to allow multiple lateral bases to be compared
        other = tools.common.convert2nparray(other)
        result = np.ones_like(other, dtype=bool)

        # check if the instances have all required metadata
        self.check()
        for indx,oth in enumerate(other):
            try:
                oth.check()
            except AttributeError: # if either is not this class instance
                result[indx] = False

            # check type
            if self.type != oth.type: result[indx] = False

            # check all keys
            for key in self.metadata.keys():
                if not np.array_equal(self.metadata[key],oth.metadata[key]): result[indx] = False
        # assume equal otherwise
        return result[0] if len(result)==1 else result

    #########################       methods       #############################

    def check(self):
        """
        Checks that object contains all attributes required for evaluating a
        particular basis set.
        """
        if self.type == 'SPHERICAL SPLINES':
            for key in ['xlaspl','xlospl','xraspl']:
                try:
                    knots = self.metadata[key]
                except KeyError:
                    print('Current attributes : ',self.metadata.keys())
                    raise KeyError('Attribute '+key+' missing for lateral basis type '+self.type)
        elif self.type == 'SPHERICAL HARMONICS':
            for key in ['lmaxhor']:
                try:
                    knots = self.metadata[key]
                except KeyError:
                    print('Current attributes : ',self.metadata.keys())
                    raise KeyError('Attribute '+key+' missing for lateral basis type '+self.type)
        elif self.type == 'PIXELS':
            for key in ['xlapix','xlopix','xsipix']:
                try:
                    knots = self.metadata[key]
                except KeyError:
                    print('Current attributes : ',self.metadata.keys())
                    raise KeyError('Attribute '+key+' missing for lateral basis type '+self.type)
        else:
            raise NotImplementedError(self.type+' has not been implemented in lateral_basis.')
        return

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
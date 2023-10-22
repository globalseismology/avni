#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard AVNI format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

import os
import numpy as np #for numerical analysis

####################### IMPORT AVNI LIBRARIES  #######################################
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
        self._name = name
        self._type = types
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def __eq__(self, other):

        if not isinstance(other,Lateral_basis): return False
        # convert to array to allow multiple lateral bases to be compared
        result = np.ones_like(other, dtype=bool)

        # check if the instances have all required metadata
        self.check()
        try:
            other.check()
        except AttributeError: # if either is not this class instance
            return False

        # check type
        if self._type != other._type: return False

        # check all keys
        for key in self.keys:
            if not np.array_equal(self[key],other.metadata[key]): return False

        # assume equal otherwise
        return True

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
    def keys(self):
        return self.metadata.keys()

    #########################       methods       #############################

    def check(self):
        """
        Checks that object contains all attributes required for evaluating a
        particular basis set.
        """
        if self._type == 'SPHERICAL SPLINES':
            for key in ['xlaspl','xlospl','xraspl']:
                try:
                    knots = self[key]
                except KeyError:
                    print('Current attributes : ',self.keys)
                    raise KeyError('Attribute '+key+' missing for lateral basis type '+self._type)
        elif self._type == 'SPHERICAL HARMONICS':
            for key in ['lmaxhor']:
                try:
                    knots = self[key]
                except KeyError:
                    print('Current attributes : ',self.keys)
                    raise KeyError('Attribute '+key+' missing for lateral basis type '+self._type)
        elif self._type == 'PIXELS':
            for key in ['xlapix','xlopix','xsipix']:
                try:
                    knots = self[key]
                except KeyError:
                    print('Current attributes : ',self.keys)
                    raise KeyError('Attribute '+key+' missing for lateral basis type '+self._type)
        else:
            raise NotImplementedError(self._type+' has not been implemented in lateral_basis.')
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
            if types[ii] not in self.proj.keys() and types[ii] != self['type']:
                self.proj[types[ii]] = {}
                self.proj[types[ii]][to_name[ii]] = {'data':None,'attributes':{}}

    def eval_lateral(self,lat,lon,store=False):
        """
        Evaluate radial basis at a depth interpolated from existing projection matrices.
        """
        if self._type == 'SPHERICAL SPLINES':
            horcof = tools.eval_splcon(lat,lon,self['xlaspl'],self['xlospl'],self['xraspl'])
        elif self._type == 'SPHERICAL HARMONICS':
            horcof = tools.eval_ylm(lat,lon, self['lmaxhor'])
        elif self._type == 'PIXELS':
            horcof = tools.eval_pixel(lat,lon, self['xlapix'],self['xlopix'],self['xsipix'])
        else:
            raise NotImplementedError(self._type+' has not been implemented in lateral_basis.eval_lateral')

        if store:
            self.data['latitude'] = lat
            self.data['longitude'] = lon
            self.data['horcof'] = horcof
        else:
            return horcof
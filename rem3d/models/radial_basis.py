#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets in the standard REM3D format."""
#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

import numpy as np #for numerical analysis

####################### IMPORT REM3D LIBRARIES  #######################################
from .. import tools
from .. import constants
#######################################################################################

# Vertical basis parameter class that defines an unique combination of functions, their radial parameterization and any scaling

class Radial_basis(object):
    '''
    A class for radial bases that defines a unique combination of parameters,
    their radial parameterization and any scaling that is used.

    Structure:
    ---------

    object.data: Contains the following fields that describe the radial basis.
        depths_in_km: depth array in km
        vercof: value of the bases evaluated at those depths
        dvercof: gradient of the bases evaluated at those depths

    object.metadata: Contains metadata for various calculations.
        name: to store a name for the radial_basis
        type: type of radial basis e.g. vbspl
        attributes: a dictionary containing variables used to define this particular type
                    e.g. knots for vbspl. Checked that these are defined using self.check.
    '''
    def __init__(self,name,types,metadata=None):
        self.data = {}
        self.data['depths_in_km'] = None
        self.data['vercof'] = None
        self.data['dvercof'] = None
        self.name = name
        self.type = types
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata
        # Check if all required atributes are available
        self.check()

    def add_attribute(self,key,value):
        """
        Add attributes needed by the radial basis

        Input parameters:
        ----------------

        key: string key name

        value: values corresponding to the key
        """
        self.metadata[key] = value

    def check(self):
        """
        Checks that object contains all attributes required for evaluating a
        particular basis set.
        """
        if self.type in ['vbspl','variable splines']:
            for key in ['knots']:
                try:
                    knots = self.metadata[key]
                except KeyError:
                    print('Current attributes : ',self.metadata.keys())
                    raise KeyError('Attribute '+key+' missing for radial basis type '+self.type)
        elif self.type in ['delta','dirac delta']:
            for key in ['info']:
                try:
                    knots = self.metadata[key]
                except KeyError:
                    print('Current attributes : ',self.metadata.keys())
                    raise KeyError('Attribute '+key+' missing for radial basis type '+self.type)
        elif self.type in ['boxcar']:
            for key in ['depthtop','depthbottom']:
                try:
                    knots = self.metadata[key]
                except KeyError:
                    print('Current attributes : ',self.metadata.keys())
                    raise KeyError('Attribute '+key+' missing for radial basis type '+self.type)
        else:
            raise TypeError('metadata type note defined in eval_radial %s' % self.type)
        return knots

    def readprojfile(self,projfile):
        """
        Reads a projection matrix file that goes between the radial bases.
        """

    def eval_radial(self,depths_in_km,store=False):
        """
        Evaluates the radial bases at various depths.

        Input parameters:
        ----------------

        depths_in_km: depths where the radial parameteriation needs to be evaluated.

        store: store them in data
        """

        # convert to numpy arrays
        depths = tools.convert2nparray(depths_in_km)

        # compute the radial parameteriation in specific depths
        if self.type in ['vbspl','variable splines']:
            knots = self.metadata['knots']
            vercof, dvercof = tools.eval_vbspl(depths,knots)
        elif self.type in ['delta','dirac delta']:
            vercof = np.ones(len(depths))
            dvercof = np.zeros(len(depths))
        elif self.type in ['boxcar','constant']:
            # convert to numpy arrays
            depthtop = tools.convert2nparray(self.metadata['depthtop'])
            depthbottom = tools.convert2nparray(self.metadata['depthbottom'])

            rtop = constants.R.to('km').magnitude - depthtop
            rbottom = constants.R.to('km').magnitude - depthbottom
            rquery = constants.R.to('km').magnitude - depths
            rrange = np.vstack((rbottom,rtop)).T
            vercof, dvercof = tools.eval_polynomial(rquery,rrange,constants.R.to('km').magnitude,types = ['CONSTANT'])
        else:
            raise TypeError('metadata type not defined in eval_radial %s' % self.type)

        # Store if needed
        if store:
            self.data['vercof'] = vercof
            self.data['dvercof'] = dvercof
            self.data['depths_in_km'] = depths
        else:
            return vercof,dvercof

    def project_boxdepth(self,depth_range):
        """
        Project from current vertical basis to a vertical boxcar basis
        depth_range is named numpy array of top and bottom depths
        """



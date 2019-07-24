#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float

import numpy as np #for numerical analysis
from scipy import sparse
import pdb

####################### IMPORT REM3D LIBRARIES  #######################################
from .lateral_basis import Lateral_basis
from .radial_basis import Radial_basis
from .common import radial_attributes
from .. import tools
#######################################################################################

# kernel set
class Kernel_set(object):
    '''
    A class for kernel sets that define the G matrix for relating data d to model m, d=Gm
    '''
    def __init__(self,dictionary):
        self.metadata ={}
        self.data = {}
        self.name = dictionary['kerstr']
        self.initialize(dictionary)
        self.extract_lateral(dictionary)
        self.extract_radial(dictionary)

    def initialize(self,dictionary,required=None,optional=None):

        #defaults
        if required is None: required = ['nmodkern','ivarkern','desckern','ncoefhor','ncoefcum','nhorpar','ihorpar','ityphpar','typehpar','numvar','varstr']
        if optional is None: optional = ['forward_modeling','scaling']

        for var in required:
            try:
                self.metadata[var] = dictionary[var]
            except:
                raise KeyError('required field '+var+' not found for kernel_set')
        for var in optional:
            try:
                self.metadata[var] = dictionary[var]
            except:
                self.metadata[var] = None

    def extract_lateral(self,dictionary):
        lateral=[]
        for ihor in np.arange(self.metadata['nhorpar']):
            types = self.metadata['typehpar'][ihor]
            metadata = {}
            metadata['ncoefhor']=dictionary['ncoefhor'][ihor]
            if 'SPHERICAL HARMONICS' in types:
                metadata['lmaxhor'] = dictionary['lmaxhor'][ihor]
            elif 'PIXELS' in types:
                for field in ['xsipix','xlapix','xlopix']:
                    metadata[field] = np.array(dictionary[field][ihor], order = 'F')
            elif 'SPHERICAL SPLINES' in types:
                for field in ['ixlspl','xlaspl','xlospl','xraspl']:
                    metadata[field] = np.array(dictionary[field][ihor], order = 'F')
            else:
                raise NotImplementedError(types+' has not been implemented in kernel_set.extract_lateral')
            lateral.append(Lateral_basis(name='HPAR'+str(ihor+1), types = types, metadata=metadata))
        self.data['lateral_basis']=np.array(lateral)

    def extract_radial(self,dictionary):
        radial={}
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        for variable in dictionary['varstr']: #loop over all variables, grabbing
            radial[variable]=[]

            ivarfind =np.where(self.metadata['varstr']==variable)[0]
            if not len(ivarfind) == 1: raise AssertionError('only one parameter can be selected in eval_kernel_set')
            findrad = np.array([(ii, dictionary['desckern'][ii]) for ii in np.arange(len(dictionary['ivarkern'])) if ivarfind[0]+1 == self.metadata['ivarkern'][ii]],dtype=dt)

            metadata = {};found = False
            types = np.unique([findrad['kernel'][ii].split(',')[-2].strip() for ii in np.arange(len(findrad))])
            if not len(types) == 1: raise AssertionError('only one types is allowed')

            for jj in np.arange(len(findrad)):
                radker = findrad['kernel'][jj]
                metadata = {};
                found = False
                if 'variable splines' in radker or 'vbspl' in radker:
                    found = True
                    metadata['knots'] = [float(findrad['kernel'][ii].split(',')[-1].split('km')[0]) for ii in np.arange(len(findrad))]
                    radial[variable].append(Radial_basis(name=radker, types = 'variable splines', metadata=metadata))

                elif 'delta' in radker or 'dirac delta' in radker:
                    found = True
                    metadata['info'] = radker.split(',')[-1]
                    radial[variable].append(Radial_basis(name=radker, types = 'dirac delta', metadata=metadata))
                elif 'boxcar' in radker or 'constant' in radker:
                    found = True
                    metadata['depthtop'] = [float(findrad['kernel'][ii].split(',')[-1].split('-')[0]) for ii in np.arange(len(findrad))]
                    metadata['depthbottom'] = [float(findrad['kernel'][ii].split(',')[-1].split('-')[1].split('km')[0]) for ii in np.arange(len(findrad))]
                    radial[variable].append(Radial_basis(name=radker, types = 'boxcar', metadata=metadata))
                if not found: raise ValueError('information not found for '+radker)
        self.data['radial_basis']=radial

    def find_radial(self,parameter):
        """
        find radial indices of a given physical parameter
        """
        # all the parameter options
        stringsum = self.metadata['varstr'][0]
        for stri in self.metadata['varstr'][1:]: stringsum= stringsum+', '+stri

        # select the radial kernels for this parameter
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        ivarfind =np.where(self.metadata['varstr']==parameter)[0]
        if not len(ivarfind) == 1: raise AssertionError('only one parameter can be selected in evaluate_bases among: '+stringsum+'. Only '+str(len(ivarfind))+' found for parameter '+parameter)
        findrad = np.array([(ii, self.metadata['desckern'][ii]) for ii in np.arange(len(self.metadata['ivarkern'])) if ivarfind[0]+1 == self.metadata['ivarkern'][ii]],dtype=dt)
        return findrad

    def evaluate_bases(self,parameter,latitude = None,longitude= None,depth_in_km = None):
        """
        depth_in_km : depth in km where the projection matrix is needed.
                      If None, returns the projection matrix for the lat/lon
                      and radial basis as a dirac delta.
        """

        # find radial indices of a given physical parameter
        findrad = self.find_radial(parameter)

        # select corresponding lateral bases
        lateral_basis = self.data['lateral_basis']
        try:
            lateral_select = lateral_basis[self.metadata['ihorpar']-1][findrad['index']]
        except:
            raise ValueError('ihorpar needs to be defined for a kernel set. The HPAR for each radial kernel')

        # check if lateral paramterization is same across all radial kernels
        if not np.all(lateral_select == lateral_select[0]): return NotImplementedError('All lateral parameterizations for a given physical parameter need to be the same. The alternative has not been implemented yet')

        # evaluate the horizontal param at these locations
        if latitude is None or longitude is None:
            horcof = None
        else:
            latitude = tools.common.convert2nparray(latitude)
            longitude = tools.common.convert2nparray(longitude)
            horcof = lateral_select[0].eval_lateral(latitude,longitude)

        #make sure only one variable is selected based on parameter input
        variables = np.unique(self.metadata['varstr'][self.metadata['ivarkern']-1][findrad['index']])
        if not len(variables) == 1: raise AssertionError('only one parameter, not '+str(len(variables))+', can be selected in evaluate_bases')        # select radial bases for this variable
        radial_select = self.data['radial_basis'][variables[0]]

        # check if lateral paramterization is same across all radial kernels
        if not np.all(radial_select == radial_select[0]): return NotImplementedError('All radial parameterizations for a given physical parameter need to be the same. The alternative has not been implemented yet')

        # evaluate the radial param at these depths
        if depth_in_km is None:
            vercof = None
        else:
            vercof, _ = radial_select[0].eval_radial(depth_in_km)
        return horcof, vercof
#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard AVNI format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float

import numpy as np #for numerical analysis
from scipy import sparse
import warnings
from configobj import ConfigObj
import os

####################### IMPORT AVNI LIBRARIES  #######################################
from .lateral_basis import Lateral_basis
from .radial_basis import Radial_basis
from .. import tools
from .. import constants
#######################################################################################

# kernel set
class Kernel_set(object):
    '''
    A class for kernel sets that define the G matrix for relating data d to model m, d=Gm
    '''
    def __init__(self,dictionary):
        self.metadata ={}
        self.data = {}
        self._name = dictionary['kerstr']
        self.get_attributes(dictionary)
        self.initialize(dictionary)
        self.extract_lateral(dictionary)
        self.extract_radial(dictionary)

    def __str__(self):
        if self._name is not None:
            output = "%s is a kernel set" % (self._name)
        else:
            output = "No kernel set has been initialized yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._name})'.format(self=self)

    def __getitem__(self,key):
        """returns metadata from key"""
        return self.metadata[key]

    def __setitem__(self,key,data):
        """sets data to key"""
        self.metadata[key] = data

    def initialize(self,dictionary,required=None,optional=None):

        #defaults
        if required is None: required = ['nmodkern','ivarkern','desckern','ncoefhor','ncoefcum','nhorpar','ihorpar','ityphpar','typehpar','numvar','varstr']
        if optional is None: optional = ['forward_modeling','scaling']

        for var in required:
            try:
                self[var] = dictionary[var]
            except:
                raise KeyError('required field '+var+' not found for kernel_set')
        for var in optional:
            try:
                self[var] = dictionary[var]
            except:
                self[var] = None
        # rest of them should be copied over
        for var in dictionary.keys():
            if var not in required+optional: self[var] = dictionary[var]

    #########################       decorators       ##########################

    @property
    def name(self):
        return self._name

    @property
    def keys(self):
        return self.metadata.keys()

    @property
    def scaling(self):
        return self['scaling'] if 'scaling' in self.keys else None

    #########################       methods       #############################

    def get_attributes(self,dictionary,parserfile=tools.get_configdir()+'/'+ constants.attributes):
        """
        Reads configuration file and get the basis attributes like knot depth for each
        parameter in modelarray list of coefficients. parser is the output from
        parser = ConfigObj(config) where config is the configuration *.ini file.
        """
        if self._name == None: raise ValueError("No kernel has been read into this Kernel_Set instance yet")

        # Read configuration file to a configuration parser. This is to make this available on the fly
        if (not os.path.isfile(parserfile)): raise IOError("Configuration file ("+parserfile+") does not exist")
        parser = ConfigObj(parserfile, unrepr=True)

        # Read the kernel set from the model3d dictionary
        kerstr=dictionary['kerstr']
        if kerstr not in parser['Kernel_Set'].keys():
            warnings.warn("Warning: kernel set "+kerstr+" not described in attribute file "+parserfile)
            return
        parser_keys = parser['Kernel_Set'][kerstr].keys()

        for key in parser_keys:
            if key=='radial_type':
                radial_knots = parser['Kernel_Set'][kerstr]['radial_knots']
                # Clustering configuration
                radial_type = parser['Kernel_Set'][kerstr]['radial_type']
                numtypes = len(radial_type)
                if numtypes != dictionary['numvar']: raise ValueError('number of radial parameterization types in attributes.ini '+str(numtypes)+' should match the number of variables on file'+str(dictionary['numvar']))

                # Loop over all kernel basis
                knot_depth=[]
                desckern = dictionary['desckern']
                # counter for dirac delta functions; should match dirac_knots in attributes
                for ii,kernel in enumerate(desckern):
                    ivar = dictionary['ivarkern'][ii]-1
                    contains=True
                    for tmpstr in radial_type[ivar].split(): contains = contains and (tmpstr in kernel)
                    if contains:
                        index = int(desckern[ii].split(',')[-1]) - 1
                        depth = radial_knots[radial_type[ivar]][index]
                    else:
                        depth=radial_knots[radial_type[ivar]]
                        warnings.warn(" Assuming a single depth for the radial kernel "+desckern[ii]+". Setting to "+str(depth)+" from attributes.ini")
                        dictionary['desckern'][ii]= kernel+', '+radial_type[ivar]+', '+str(depth)+' km'
                    if str(depth) not in dictionary['desckern'][ii]:
                        dictionary['desckern'][ii]=kernel+', '+str(depth)+' km'
                    knot_depth.append(depth)
                dictionary['knots']=np.array(knot_depth)
            else:
                dictionary[key]=parser['Kernel_Set'][kerstr][key]
        return

    def extract_lateral(self,dictionary):
        lateral=[]
        for ihor in np.arange(self['nhorpar']):
            types = self['typehpar'][ihor]
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
                raise ValueError(types+' has not been implemented in kernel_set.extract_lateral')
            lateral.append(Lateral_basis(name='HPAR'+str(ihor+1), types = types, metadata=metadata))
        self.data['lateral_basis']=np.array(lateral)

    def extract_radial(self,dictionary):
        radial={}
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        for variable in dictionary['varstr']: #loop over all variables, grabbing
            radial[variable]=[]
            ivarfind =np.where(self['varstr']==variable)[0]
            if not len(ivarfind) == 1: raise AssertionError('only one parameter can be selected in eval_kernel_set')
            findrad = np.array([(ii, dictionary['desckern'][ii]) for ii in np.arange(len(dictionary['ivarkern'])) if ivarfind[0]+1 == self['ivarkern'][ii]],dtype=dt)

            metadata = {};found = False
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
                if not found:
                    raise ValueError('information not found for lateral kernel '+radker)
        self.data['radial_basis']=radial

    def find_radial(self,parameter):
        """
        find radial indices of a given physical parameter
        """
        # all the parameter options
        stringsum = self['varstr'][0]
        for stri in self['varstr'][1:]: stringsum= stringsum+', '+stri

        # select the radial kernels for this parameter
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        ivarfind =np.where(self['varstr']==parameter)[0]
        if not len(ivarfind) == 1: raise AssertionError('only one parameter can be selected in evaluate_bases among: '+stringsum+'. Only '+str(len(ivarfind))+' found for parameter '+parameter)
        findrad = np.array([(ii, self['desckern'][ii]) for ii in np.arange(len(self['ivarkern'])) if ivarfind[0]+1 == self['ivarkern'][ii]],dtype=dt)
        return findrad

    def search_radial(self,search,unique=False):
        found=False
        output=[]
        for radker in self['varstr']: #loop over all variables, grabbing
            if search.replace(" ", "").lower() in radker.replace(" ", "").lower():
                if found and unique: raise ValueError('String '+search+' maps to multiple radial kernels. Please be more specific in attributes.ini')
                found=True
                output.append(radker)
        return output[0] if unique else output

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
            lateral_select = lateral_basis[self['ihorpar']-1][findrad['index']]
        except:
            raise ValueError('ihorpar needs to be defined for a kernel set. The HPAR for each radial kernel')

        # check if lateral paramterization is same across all radial kernels
        for lateral in lateral_select:
            if not np.all(lateral == lateral_select[0]):
                return NotImplementedError('All lateral parameterizations for a given physical parameter need to be the same. The alternative has not been implemented yet')

        # evaluate the horizontal param at these locations
        if latitude is None or longitude is None:
            horcof = None
        else:
            latitude = tools.common.convert2nparray(latitude)
            longitude = tools.common.convert2nparray(longitude)
            horcof = lateral_select[0].eval_lateral(latitude,longitude)

        #make sure only one variable is selected based on parameter input
        variables = np.unique(self['varstr'][self['ivarkern']-1][findrad['index']])
        if not len(variables) == 1: raise AssertionError('only one parameter, not '+str(len(variables))+', can be selected in evaluate_bases')        # select radial bases for this variable
        radial_select = self.data['radial_basis'][variables[0]]

        # check if lateral paramterization is same across all radial kernels
        for radial in radial_select:
            if not np.all(radial == radial_select[0]):
                return NotImplementedError('All radial parameterizations for a given physical parameter need to be the same. The alternative has not been implemented yet')
            if 'delta' not in radial._type and depth_in_km is None:
                raise ValueError('depth_in_km needs to be specified for evaluate_bases if parameterization is not dirac delta  i.e. '+radial._type)

        # evaluate the radial param at these depths
        if depth_in_km is None:
            vercof = None
        else:
            vercof, _ = radial_select[0].eval_radial(depth_in_km)
        return horcof, vercof

    def pixeldepths(self,parameter):
        typehpar = self.data['lateral_basis']
        if not len(typehpar) == 1: raise AssertionError('only one type of horizontal parameterization allowed',)
        if not typehpar[0].type == 'PIXELS': raise AssertionError('Only PIXELS allowed, not '+typehpar[0].type)
        kernel_param = self.data['radial_basis'][parameter]
        depths = []
        for index,radker in enumerate(kernel_param):
            if  'depthtop' in radker.keys and 'depthbottom' in radker.keys:
                depths.append((radker.metadata['depthtop'][index]+radker.metadata['depthbottom'][index])/2.)
        return np.asarray(depths) if len(depths) > 0 else None
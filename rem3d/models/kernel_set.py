#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import numpy as np #for numerical analysis
from datetime import date  #to give a timestamp to output and compare times
import pdb    #for the debugger pdb.set_trace()

####################### IMPORT REM3D LIBRARIES  #######################################
from .lateral_basis import lateral_basis
from .radial_basis import radial_basis
from .common import radial_attributes
from ..f2py import splcon
#######################################################################################

# kernel set
class kernel_set(object):
    '''
    A class for kernel sets that define the G matrix for relating data d to model m, d=Gm
    '''
    def __init__(self,dict):
        self.metadata ={}
        self.data = {}
        self.name = dict['kerstr']
        self.initialize(dict)     
        self.extract_lateral(dict)
        #self.extract_radial(dict)
#         names, types, attributes = radial_attributes(dict['desckern'])
#         for name in names:
#             radial_basis(name,type,attributes)
#         
#         #dict['varstr']
#         # intialize lateral basis
#         for ihor in range(dict['nhorpar']): 
#             lateral_basis(dict['typehpar'][ihor])
            
    def initialize(self,dict,required = ['nmodkern','ivarkern','desckern','ncoefhor','ncoefcum','nhorpar','ihorpar','ityphpar','typehpar','numvar','varstr'],optional = ['forward_modeling','scaling']):
        for var in required:
            try:
                self.metadata[var] = dict[var]
            except:
                raise KeyError('required field '+var+' not found for kernel_set')
        for var in optional:
            try:
                self.metadata[var] = dict[var]
            except:
                self.metadata[var] = None
        
    def extract_lateral(self,dict):
        lateral=[]
        for ihor in range(self.metadata['nhorpar']):
            type = self.metadata['typehpar'][ihor]
            attributes = {}
            attributes['ncoefhor']=dict['ncoefhor'][ihor]
            if 'SPHERICAL HARMONICS' in type:
                attributes['lmaxhor'] = dict['lmaxhor'][ihor]
            elif 'PIXELS' in type:
                for field in ['xsipix','xlapix','xlopix']:
                    attributes[field] = np.array(dict[field][ihor], order = 'F')
            elif 'SPHERICAL SPLINES' in type:
                for field in ['ixlspl','xlaspl','xlospl','xraspl']:
                    attributes[field] = np.array(dict[field][ihor], order = 'F')
            pdb.set_trace()
            #ncon,icon,con = splcon(0.,0.,attributes['ncoefhor'],attributes['xlaspl'],attributes['xlospl'],attributes['xraspl'])
            lateral.append(lateral_basis(name='HPAR'+str(ihor+1), type = type, attributes=attributes))

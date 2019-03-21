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
        self.metadata['scaling'] = dict['scaling']
        #self.metadata['forward_modeling'] = dict['forward_modeling']
        self.metadata['desckern'] = dict['desckern']
        pdb.set_trace()
        
        names, types, attributes = radial_attributes(dict['desckern'])
        for name in names:
            radial_basis(name,type,attributes)
        
        #dict['varstr']
        # intialize lateral basis
        for ihor in range(dict['nhorpar']): 
            lateral_basis(dict['typehpar'][ihor])
            
            

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

#######################################################################################

# kernel set
class kernel_set(object):
    '''
    A class for kernel sets that define the G matrix for relating data d to model m, d=Gm
    '''

    def __init__(self,parameters,radial_basis,lateral_basis,indices):
        self.metadata ={}
        self.data = {}
        self.name = None
        self.type = None
        self.refmodel = None
        self.description = None
        if file is not None: self.readfile(file)



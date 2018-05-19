#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format. 
Author: Raj Moulik, 2018"""

#####################  IMPORT STANDARD MODULES   ######################################   

import sys,os
import argparse #parsing arguments
import glob # pattern matching for finding files
import numpy as np #for numerical analysis
from scipy import stats #for stats analysis
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output and compare times
import pdb    #for the debugger pdb.set_trace()
from math import pi
import fortranformat as ff #reading/writing fortran formatted text
from future.utils import native_str

#####################

def readepixfile(filename):
    """Read .epix file format."""

    currentdir=os.getcwd()
    try: 
        f = open(filename, 'r')
    except IOError:
        print "File ("+filename+") does not exist in the current directory - "+currentdir
        sys.exit(2)
            
    epixarr=np.genfromtxt(filename, dtype=None,comments="#",names=['lat','lon','pixsize','val'])

    return epixarr

Model1D_Attr = np.dtype([
    (native_str('radius'),np.float),
    (native_str('rho'),np.float),
    (native_str('vp'),np.float),
    (native_str('vph'),np.float),
    (native_str('vpv'),np.float),
    (native_str('vs'),np.float),
    (native_str('vsh'),np.float),
    (native_str('vsv'),np.float),
    (native_str('Qmu'),np.float),
    (native_str('Qkappa'),np.float),
    (native_str('eta'),np.float)])

class reference1D(object):
    '''
    A class for 1D reference Earth models used in tomography
    '''

    def __init__(self):

        self.__nlayers__ = None
	self.data = None
	self.metdata = None

    def read(self,path_to_model,fmt='card'):

        if fmt=='card':
            modelarr = np.genfromtxt(path_to_model,dtype=None,comments='#',skip_header=3,
	        names=['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta'])
      
        else:
	    raise ValueError('model format ',fmt,' is not currently implemented')

	self.__nlayers__ = len(modelarr['radius'])

	#Do voigt averaging here?
	vs = np.zeros(len(modelarr['radius']))
	vp = np.zeros(len(modelarr['radius']))

	self.data = np.empty(self.__nlayers__,dtype=Model1D_Attr)
	self.data['radius'] = modelarr['radius']
	self.data['rho'] = modelarr['rho']
	self.data['vpv'] = modelarr['vpv']
	self.data['vsv'] = modelarr['vsv']
	self.data['Qkappa'] = modelarr['Qkappa']
	self.data['Qmu'] = modelarr['Qmu']
	self.data['vph'] = modelarr['vph']
	self.data['vsh'] = modelarr['vsh']
	self.data['eta'] = modelarr['eta']
	self.data['vp'] = vp
	self.data['vs'] = vs

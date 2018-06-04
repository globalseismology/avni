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
	self.name = None
	self.radius_max = None

    def read(self,path_to_model,fmt='card'):

        if fmt=='card':
            modelarr = np.genfromtxt(path_to_model,dtype=None,comments='#',skip_header=3,
	        names=['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta'])
      
        else:
	    raise ValueError('model format ',fmt,' is not currently implemented')

	self.__nlayers__ = len(modelarr['radius'])

	#Do voigt averaging here?
	#vp = voigt_average(modelarr['vph'],modelarr['vpv'])
	#vs = voigt_average(modelarr['vsh'],modelarr['vsv'])
	vp = np.zeros(len(modelarr['radius']))
	vs = np.zeros(len(modelarr['radius']))

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
        self.radius_max = np.max(self.data['radius'])

    def to_TauPmodel(self,model_name,fmt='tvel'):
        '''
        Writes a model file that is compatible with TauP.
	file format options 'tvel' and 'nd'.

	Note: TauP can't handle zero shear velocity in the ocean layer...
	      To work around this, zero values an ocean layer will be written 
	      as 1e-4.
        '''

	f = open(model_name+'.tvel','w')
	f.write('{} - P\n'.format(model_name))
	f.write('{} - S\n'.format(model_name))

        for i in range(0,len(self.data)):
	    f.write('{:2.4f}   {:2.4f}   {:2.4f}    {:2.4f}\n'.format(
	    (self.radius_max - self.data['radius'][::-1][i]) / 1000.0,
	    self.data['vp'][::-1][i] / 1000.0,
	    self.data['vs'][::-1][i] / 1000.0,
            self.data['rho'][::-1][i] / 1000.0))
	    
	f.close()

    def to_axisem(self,model_name,anelastic=True,anisotropic=True):
        '''
	Write 1D model to be used as an external model in axisem
        '''
	f = open(model_name+'.bm','w')
	n_discon = 1

	if anelastic:
	    f.write('ANELASTIC     T\n')
	else:
	    f.write('ANELASTIC     F\n')

	if anisotropic:
	    f.write('ANISOTROPIC     T\n')
	else:
	    f.write('ANISOTROPIC     F\n')

	f.write('UNITS      m\n')

	if anisotropic:
	   f.write('COLUMNS   radius    rho    vpv    vsv    qka    qmu    vph    vsh    eta\n')

	   for i in range(0,len(self.data)):
	       f.write('{}    {}    {}    {}    {}    {}    {}    {}    {}\n'.format(
	       self.data['radius'][::-1][i],
	       self.data['rho'][::-1][i],
	       self.data['vpv'][::-1][i],
	       self.data['vsv'][::-1][i],
	       self.data['Qkappa'][::-1][i],
	       self.data['Qmu'][::-1][i],
	       self.data['vph'][::-1][i],
	       self.data['vsh'][::-1][i],
	       self.data['eta'][::-1][i]) )

	       if i < len(self.data)-1 and self.data['radius'][::-1][i] == self.data['radius'][::-1][i+1]:
	           depth_here = (self.radius_max - self.data['radius'][::-1][i]) / 1000.0 
	           f.write('#    Discontinuity {}, depth {:6.2f} km\n'.format(n_discon,depth_here))
		   n_discon += 1

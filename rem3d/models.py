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
    """Read .epix file format from a file.

    Parameters
    ----------

    filename : Name of the file containing four columns
              (latitude, longitude, pixel_size, value)

    """

    currentdir=os.getcwd()
    try:
        f = open(filename, 'r')
        epixarr=np.genfromtxt(filename, dtype=None,comments="#",names=['lat','lon','pixsize','val'])
        ierror=0
    except IOError:
        sys.exit("File ("+filename+") does not exist in the current directory - "+currentdir)

    return ierror,epixarr

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
    (native_str('eta'),np.float),
    (native_str('A'),np.float),
    (native_str('C'),np.float),
    (native_str('N'),np.float),
    (native_str('L'),np.float),
    (native_str('F'),np.float),
    (native_str('kappa'),np.float),
    (native_str('mu'),np.float)])

class reference1D(object):
    '''
    A class for 1D reference Earth models used in tomography
    '''

    def __init__(self):
        self.__nlayers__ = None
        self.data = None
        self.metadata = None
        self.name = None
        self.radius_max = None

    def read(self,file,fmt='card'):
        '''
        Read a card deck file used in OBANI. Other formats not ready yet
        '''
        if fmt=='card':
            modelarr = np.genfromtxt(file,dtype=None,comments='#',skip_header=3,
            names=['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta'])
        else:
            raise ValueError('model format ',fmt,' is not currently implemented')

        self.__nlayers__ = len(modelarr['radius'])

        self.data = np.zeros(self.__nlayers__,dtype=Model1D_Attr)
        self.data['radius'] = modelarr['radius']
        self.data['rho'] = modelarr['rho']
        self.data['vpv'] = modelarr['vpv']
        self.data['vsv'] = modelarr['vsv']
        self.data['Qkappa'] = modelarr['Qkappa']
        self.data['Qmu'] = modelarr['Qmu']
        self.data['vph'] = modelarr['vph']
        self.data['vsh'] = modelarr['vsh']
        self.data['eta'] = modelarr['eta']
        self.radius_max = np.max(self.data['radius'])

    def get_Love_elastic(self):
        '''
        Get the Love parameters and Voigt averaged elastic properties with depth
        '''
        if self.data is not None and self.__nlayers__ > 0:
            self.data['A'] = self.data['rho']*self.data['vph']**2
            self.data['C'] = self.data['rho']*self.data['vpv']**2
            self.data['N'] = self.data['rho']*self.data['vsh']**2 
            self.data['L'] = self.data['rho']*self.data['vsv']**2
            self.data['F'] = self.data['eta']*(self.data['A']-2.*self.data['L'])
            self.data['kappa'] = (4.0*(self.data['A']+self.data['F']-self.data['N'])+self.data['C'])/9.
            self.data['mu'] = (self.data['A']+self.data['C']-2.*self.data['F']+5.*self.data['N']+6.*self.data['L'])/15.
            self.data['vp'] = np.sqrt(np.divide((self.data['kappa']+4.*self.data['mu']/3.),self.data['rho']))
            self.data['vs'] = np.sqrt(np.divide(self.data['mu'],self.data['rho']))

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
        n_discon = 0

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
                n_discon += 1
                f.write('#    Discontinuity {}, depth {:6.2f} km\n'.format(n_discon,depth_here))

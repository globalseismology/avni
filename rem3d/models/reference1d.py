#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int,list

import numpy as np #for numerical analysis
import fortranformat as ff #reading/writing fortran formatted text
from six import string_types # to check if variable is string using isinstance
from numpy.lib.recfunctions import append_fields # for appending named columns to named numpy arrays
from scipy.interpolate import griddata
from copy import deepcopy
from collections import Counter
import traceback
import pandas as pd
import pdb
import pint
import ntpath

####################### IMPORT REM3D LIBRARIES  #######################################
from .. import constants
from .. import tools
from rem3d.f2py import getbullen
#######################################################################################
# 1D model class

class Reference1D(object):
    '''
    A class for 1D reference Earth models used in tomography
    '''

    def __init__(self,file=None):
        self.data = None
        self.metadata = {}
        self.name = None
        self._radius_max = None
        self._nlayers = None
        if file is not None:
            self.read(file)
            self.get_Love_elastic()
            self.get_discontinuity()
            self.get_mineralogical()

    def __str__(self):
        if self.data is not None and self._nlayers > 0:
            output = "%s is a one-dimensional model with %s layers and radius up to %s km" % (self.name, self._nlayers,self._radius_max/1000.)
        else:
            output = "No model has been read into this reference1D instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self.name,self._radius_max})'.format(self=self)

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result

    def read(self,file):
        '''
        Read a card deck file used in OBANI. Other formats not ready yet
        '''
        try:
            self.readmineoscards(file)
        except:
            pdb.set_trace()
            var1 = traceback.format_exc()
            print(var1)
            raise NotImplementedError('model format is not currently implemented in reference1D.read')

    def readmineoscards(self,file):
        # Operations between PintArrays of different unit registry will not work.
        # We can change the unit registry that will be used in creating new
        # PintArrays to prevent this issue.
        pint.PintType.ureg = constants.ureg

        names=['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta']
        units =['m','kg/m^3','m/s','m/s','dimensionless','dimensionless','m/s','m/s','dimensionless']
        fields=list(zip(names,units))
        #formats=[np.float for ii in range(len(fields))]
        # modelarr = np.genfromtxt(file,dtype=None,comments='#',skip_header=3,names=fields)
        modelarr = pd.read_csv(file,skiprows=3,comment='#',sep='\s+',names=fields)
        # read the punit units from last header
        modelarr_ = modelarr.pint.quantify(level=-1)
        self.metadata['attributes'] = names
        self.metadata['description'] = 'Read from '+file
        self.metadata['filename'] = file
        self.name = ntpath.basename(file)
        self._nlayers = len(modelarr['radius'])
        # Create data array
        PA_ = pint.PintArray
        modelarr_['depth'] = PA_((constants.R.magnitude - modelarr_['radius'].pint.to(constants.R.units).data).tolist(), dtype = constants.R.units)
        self.data = modelarr_
        self._radius_max = max(self.data['radius']).magnitude

    def get_Love_elastic(self):
        '''
        Get the Love parameters and Voigt averaged elastic properties with depth

        A,C,N,L,F: anisotropy elastic Love parameters

        kappa: bulk modulus

        mu: shear modulus

        vphi: bulk sound velocity

        xi: shear anisotropy ratio

        phi: P anisotropy ratio

        Zs, Zp: S and P impedances
        '''
        if self.data is not None and self._nlayers > 0:
            # Add metadata
            for field in ['A','C','N','L','F','vp','vs','vphi','xi','phi','Zp','Zs']: self.metadata['attributes'].append(field)

            # Add data fields
            self.data['A'] = self.data['rho']*self.data['vph']**2
            self.data['C'] = self.data['rho']*self.data['vpv']**2
            self.data['N'] = self.data['rho']*self.data['vsh']**2
            self.data['L'] = self.data['rho']*self.data['vsv']**2
            self.data['F'] = self.data['eta']*(self.data['A']-2.*self.data['L'])

            # equivalent isotropic
            self.data['kappa'] = (4.0*(self.data['A']+self.data['F']-self.data['N'])+self.data['C'])/9.
            self.data['mu'] = (self.data['A']+self.data['C']-2.*self.data['F']+5.*self.data['N']+6.*self.data['L'])/15.
            self.data['vp'] = ((self.data['kappa']+4.*self.data['mu']/3.)/self.data['rho']).pow(0.5)
            self.data['vs'] = (self.data['mu']/self.data['rho']).pow(0.5)
            self.data['vphi'] = (self.data['kappa']/self.data['rho']).pow(0.5)

            # anisotropy
            self.data['xi'] = (self.data['vsh'].div(self.data['vsv'])).pow(2)
            self.data['phi'] = (self.data['vpv'].div(self.data['vph'])).pow(2)

            # impedance contrasts
            self.data['Zp'] = self.data['vp']*self.data['rho']
            self.data['Zs'] = self.data['vs']*self.data['rho']
        else:
            raise ValueError('reference1D object is not allocated')

    def get_mineralogical(self):
        '''
        Get the Love parameters and Voigt averaged elastic properties with depth

        gravity: gavity at each depth

        Brunt-Vaisala Frequency: Used for Bullen's parameter

        Bullen: Bullen's parameter

        pressure: pressure at each depth
        '''
        if self.data is not None and self._nlayers > 0:
            if constants.planetpreferred == 'Earth':
                file = self.metadata['filename']
                layers = self._nlayers
                grav,vaisala,bullen,pressure = getbullen(file,layers,constants.omega.to_base_units().magnitude,constants.G.to_base_units().magnitude)
                # Add metadata
                for field in ['gravity','Brunt-Vaisala','Bullen','pressure']: self.metadata['attributes'].append(field)

                # Add data fields
                pdb.set_trace()
                self.data=append_fields(self.data, 'gravity', grav, usemask=False)
                self.data=append_fields(self.data, 'Brunt-Vaisala', vaisala, usemask=False)
                self.data=append_fields(self.data, 'Bullen', bullen, usemask=False)
                self.data=append_fields(self.data, 'pressure', pressure, usemask=False)
            else:
                print('Warning: mineralogical parameters not evaluated for '+constants.planetpreferred)
        else:
            raise ValueError('reference1D object is not allocated')

    def get_discontinuity(self):
        '''
        Get values, average values and contrasts at discontinuities

        Output:
        ------

        Returns a structure self.metadata['disc'] that has three arrays:

        delta: containing absolute difference in attributes between smaller/larger radii

        average: containing absolute average attributes between smaller/larger radii

        contrasts: containing contrast in attributes (in %)
        '''
        disc_depths = [item.magnitude for item, count in Counter(self.data['depth']).items() if count > 1]
        disc = {}
# Create a named array for discontinuities

        for field in ['delta','average','contrast']: disc[field] = self.data.copy().drop(range(len(np.unique(disc_depths)),len(self.data)))

        # default names and units as percent
        names = self.data.columns.tolist()
        units = ['percent' for name in names]
        fields=list(zip(names,units))

        for icount,depth in enumerate(disc_depths):
            sel = self.data[self.data['depth'].data==depth]
            for field in sel:
                if field == 'radius' or field == 'depth':
                    disc['delta'][field][icount] = sel[field].iat[0]
                    disc['average'][field][icount] = sel[field].iat[0]
                    pdb.set_trace()
                    disc['contrast'][field][icount] = sel[field].iat[0]
                else:
                    disc['delta'][field][icount] = sel[field].iat[0]-sel[field].iat[1]
                    disc['average'][field][icount] = 0.5*(sel[field].iat[0]+sel[field].iat[1])
                    pdb.set_trace()
                    ## contrasts need to be in %
                    contrast = (abs(disc['delta'][field][icount]) / disc['average'][field][icount]).to('percent')
                    disc['contrast'][field][icount] = (abs(disc['delta'][field][icount]) / disc['average'][field][icount]).to('percent')


        #---- try to find discontinuities
        discfind = disc['delta']['radius'][np.abs(1221.5-disc['delta']['radius']/1000.)<25.]
        if len(discfind) <= 0: # not found
            print("Warning: itopic not found")
        elif len(discfind) > 1: raise ValueError('get_discontinuity: multiple values within discontinuity limits')
        else:
            disc['itopic'] = np.where(self.data['radius']==discfind[0])[0][1]

        discfind = disc['delta']['radius'][np.abs(3480.0-disc['delta']['radius']/1000.)<25.]
        if len(discfind) <= 0: # not found
            print("Warning: itopoc not found")
        elif len(discfind) > 1:
            raise ValueError('get_discontinuity: multiple values within discontinuity limits')
        else:
            disc['itopoc'] = np.where(self.data['radius']==discfind[0])[0][1]

        ###   Top of crust
        discfind = np.where(np.logical_and(self.data['vp']<7500.,self.data['vs']>0.))[0]
        if len(discfind) > 0: disc['itopcrust'] = max(discfind) + 1
        #discfind = disc['delta']['radius'][np.abs(6368.0-disc['delta']['radius']/1000.)<0.1]
#         if len(discfind) <= 0: # not found
#             print("Warning: itopcrust not found")
#         elif len(discfind) > 1:
#             raise ValueError('get_discontinuity: multiple values within discontinuity limits')
#         else:
            #disc['itopcrust'] = np.where(self.data['radius']==discfind[0])[0][1]

        itopmantle = min(np.where(self.data['vp']<7500.)[0])
        if itopmantle >0: disc['itopmantle'] = itopmantle

        self.metadata['discontinuities'] = disc

    def get_custom_parameter(self,parameters):
        '''
        Get the arrays of custom parameters defined in various Earth models
        '''
        if self.data is not None and self._nlayers > 0:
            # convert to array for ease of looping
            if isinstance(parameters,string_types): parameters = np.array([parameters])

            for ii in np.arange(parameters.size):
                if parameters[ii] not in list(self.data.dtype.names):
                    if 'as' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], np.divide(self.data['vsh'] - self.data['vsv'],self.data['vs'],out=np.zeros_like(self.data['vs']), where= self.data['vs'] != 0.)*100. , usemask=False)
                    elif 'ap' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], np.divide(self.data['vph'] - self.data['vpv'],self.data['vp'],out=np.zeros_like(self.data['vp']), where= self.data['vp'] != 0.)*100. , usemask=False)
                    else:
                        raise NotImplementedError('parameter ',parameters[ii],' is not currently implemented in reference1D.get_custom_parameter')
        else:
            raise ValueError('reference1D object is not allocated')

    def evaluate_at_depth(self,depth_in_km,parameter='vs',interpolation='linear'):
        '''
        Get the values of a parameter at a given depth
        '''
        values=None
        depth_in_km = tools.convert2nparray(depth_in_km)

        if self.data is not None and self._nlayers > 0:
            if parameter in self.data.dtype.names:
                values = self.data[parameter]
                depth_array = (constants.R.to_base_units().magnitude - self.data['radius'])/1000. # in km
                # Sort to make interpolation possible
                indx = depth_array.argsort()
                values = griddata(depth_array[indx], values[indx], depth_in_km, method=interpolation)
                if len(depth_in_km)==1: values = values.item()
            else:
                raise ValueError('parameter '+parameter+' not defined in array')
        else:
            raise ValueError('reference1D object is not allocated')
        return values

    def to_mineoscards(self,directory='.',fmt='cards'):
        '''
        Writes a model file that is compatible with MINEOS.
        '''
        parameters = ['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta']
        if self.data is not None and self._nlayers > 0:
            model_name = self.name
            ntotlev = self._nlayers
            itopic = self.metadata['discontinuities']['itopic']
            itopoc = self.metadata['discontinuities']['itopoc']
            itopmantle = self.metadata['discontinuities']['itopmantle']
            itopcrust = self.metadata['discontinuities']['itopcrust']

            f = open(directory+'/'+model_name+'.'+fmt,'w')
            f.write(model_name+'\n')
            f.write('1 1. 1 1\n')
            line = ff.FortranRecordWriter('(5I5)')
            f.write(line.write([ntotlev,itopic,itopoc,itopmantle,itopcrust])+u'\n')
            line = ff.FortranRecordWriter('(f8.0,3f9.2,2f9.1,2f9.2,f9.5)')

            write = self.data[parameters]
            for _,val in enumerate(write):
                f.write(line.write(val)+u'\n')
            f.close()
        else:
            raise ValueError('reference1D object is not allocated')


    def to_TauPmodel(self,directory='.',fmt='tvel'):
        '''
        Writes a model file that is compatible with TauP.
        file format options 'tvel' and 'nd'.

        Note: TauP can't handle zero shear velocity in the ocean layer...
          To work around this, zero values an ocean layer will be written
          as 1e-4.
        '''
        if self.data is not None and self._nlayers > 0:
            model_name = self.name
            f = open(directory+'/'+model_name+'.'+fmt,'w')
            f.write('{} - P\n'.format(model_name))
            f.write('{} - S\n'.format(model_name))

            for i in range(0,len(self.data)):
                f.write('{:2.4f}   {:2.4f}   {:2.4f}    {:2.4f}\n'.format(
                   (self._radius_max - self.data['radius'][::-1][i]) / 1000.0,
                   self.data['vp'][::-1][i] / 1000.0,
                   self.data['vs'][::-1][i] / 1000.0,
                   self.data['rho'][::-1][i] / 1000.0))
            f.close()
        else:
            raise ValueError('reference1D object is not allocated')

        if self.data['vp'][-1] == 0 or self.data['vs'][-1] == 0:
            raise Warning('zero velocity layer detected at surface ...\n \
                      TauP raytracing may not work')

    def to_axisem(self,directory='.',anelastic=True,anisotropic=True):
        '''
         Write 1D model to be used as an external model in axisem
        '''
        if self.data is not None and self._nlayers > 0:
            model_name = self.name
            f = open(directory+'/'+model_name+'.bm','w')
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
                    depth_here = (self._radius_max - self.data['radius'][::-1][i]) / 1000.0
                    n_discon += 1
                    f.write('#    Discontinuity {}, depth {:6.2f} km\n'.format(n_discon,depth_here))
        else:
            raise ValueError('reference1D object is not allocated')
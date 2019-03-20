#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import numpy as np #for numerical analysis
import pdb    #for the debugger pdb.set_trace()
import fortranformat as ff #reading/writing fortran formatted text
from future.utils import native_str
from six import string_types # to check if variable is string using isinstance
from numpy.lib.recfunctions import append_fields # for appending named columns to named numpy arrays
from scipy.interpolate import griddata
import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases
import matplotlib.pyplot as plt
from matplotlib import gridspec # Relative size of subplots
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
from copy import copy, deepcopy
from collections import Counter

####################### IMPORT REM3D LIBRARIES  #######################################
from .. import plots 
from .. import constants
#######################################################################################
# 1D model class

class reference1D(object):
    '''
    A class for 1D reference Earth models used in tomography
    '''

    def __init__(self,file=None):
        self.__nlayers__ = None
        self.data = None
        self.metadata = {}
        self.name = None
        self.radius_max = None
        if file is not None: 
            self.read(file)
            self.get_Love_elastic()
            self.get_discontinuity()
    
    def __str__(self):
        if self.data is not None and self.__nlayers__ > 0:
            output = "%s is a one-dimensional model with %s layers and radius up to %s km" % (self.name, self.__nlayers__,self.radius_max/1000.)
        else:
            output = "No model has been read into this reference1D instance yet"
        return output
        
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
            
    def read(self,file,fmt='card'):
        '''
        Read a card deck file used in OBANI. Other formats not ready yet
        '''
        if fmt=='card':
            names=['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta']
            formats=[np.float for ii in range(len(names))]
            modelarr = np.genfromtxt(file,dtype=None,comments='#',skip_header=3,
            names=names)
            # Add depth assuming model describes from Earth center to surface
            names.append('depth'); formats.append(np.float)
            modelarr=append_fields(modelarr, 'depth', constants.R - modelarr['radius'], usemask=False)
            self.metadata['attributes'] = names
            self.metadata['description'] = 'Read from '+file
            self.metadata['filename'] = file
            self.name = ntpath.basename(file)
        else:
            raise NotImplementedError('model format ',fmt,' is not currently implemented in reference1D.read')

        self.__nlayers__ = len(modelarr['radius'])
        # Create data array
        Model1D_Attr = np.dtype([(native_str(names[ii]),formats[ii]) for ii in range(len(names))])
        self.data = np.zeros(self.__nlayers__,dtype=Model1D_Attr)
        self.data['radius'] = modelarr['radius']
        self.data['depth'] = modelarr['depth']
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
        
        A,C,N,L,F: anisotropy elastic Love parameters
        kappa: bulk modulus
        mu: shear modulus
        vphi: bulk sound velocity
        xi: shear anisotropy ratio
        phi: P anisotropy ratio
        Zs, Zp: S and P impedances
        '''
        if self.data is not None and self.__nlayers__ > 0:
            # Add metadata
            self.metadata['attributes'].append(['A','C','N','L','F','vp','vs','vphi','xi','phi','Zp','Zs'])
            
            # Add data fields
            self.data=append_fields(self.data, 'A', self.data['rho']*self.data['vph']**2 , usemask=False)
            self.data=append_fields(self.data, 'C', self.data['rho']*self.data['vpv']**2 , usemask=False)
            self.data=append_fields(self.data, 'N', self.data['rho']*self.data['vsh']**2 , usemask=False)
            self.data=append_fields(self.data, 'L', self.data['rho']*self.data['vsv']**2 , usemask=False)
            self.data=append_fields(self.data, 'F', self.data['eta']*(self.data['A']-2.*self.data['L']) , usemask=False)
            self.data=append_fields(self.data, 'kappa', (4.0*(self.data['A']+self.data['F']-self.data['N'])+self.data['C'])/9. , usemask=False)
            self.data=append_fields(self.data, 'mu', (self.data['A']+self.data['C']-2.*self.data['F']+5.*self.data['N']+6.*self.data['L'])/15. , usemask=False)
            self.data=append_fields(self.data, 'vp', np.sqrt(np.divide((self.data['kappa']+4.*self.data['mu']/3.),self.data['rho'])) , usemask=False)
            self.data=append_fields(self.data, 'vs', np.sqrt(np.divide(self.data['mu'],self.data['rho'])) , usemask=False)
            self.data=append_fields(self.data, 'vphi', np.sqrt(np.divide(self.data['kappa'],self.data['rho'])) , usemask=False)
            with np.errstate(divide='ignore', invalid='ignore'): # Ignore warning about dividing by zero
                xi = np.power(np.divide(self.data['vsh'],self.data['vsv']),2)
            self.data=append_fields(self.data, 'xi', xi , usemask=False)
            self.data=append_fields(self.data, 'phi', np.power(np.divide(self.data['vpv'],self.data['vph']),2) , usemask=False)
            self.data=append_fields(self.data, 'Zp', self.data['vp']*self.data['rho'], usemask=False)
            self.data=append_fields(self.data, 'Zs', self.data['vs']*self.data['rho'], usemask=False)
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
        
        disc_depths = [item for item, count in Counter(self.data['depth']).items() if count > 1]
        disc = {}
# Create a named array for discontinuities
        disc['delta'] = np.zeros(len(np.unique(disc_depths)),dtype=self.data.dtype)
        disc['contrast'] = np.copy(disc['delta']);disc['average'] = np.copy(disc['delta'])
        
        icount  = 0
        for depth in np.unique(disc_depths):
            sel = self.data[np.where(self.data['depth']==depth)]
            for field in self.data.dtype.names:
                if field == 'radius' or field == 'depth':
                    disc['delta'][field][icount] = sel[0][field]
                    disc['average'][field][icount] = sel[0][field]
                    disc['contrast'][field][icount] = sel[0][field]
                else:
                    disc['delta'][field][icount] = sel[0][field]-sel[1][field] 
                    disc['average'][field][icount] = 0.5*(sel[0][field]+sel[1][field])
                    disc['contrast'][field][icount] = abs(disc['delta'][field][icount]) / disc['average'][field][icount]*100.
            icount = icount+1
            
            
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
        if self.data is not None and self.__nlayers__ > 0:
            # convert to array for ease of looping
            if isinstance(parameters,string_types): parameters = np.array([parameters]) 
            
            for ii in np.arange(parameters.size):
                if parameters[ii] in list(self.data.dtype.names):
                    print('... parameter '+parameters[ii]+' already exists in '+self.name)
                else:
                    if 'SH-SV' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], self.data['vsh'] - self.data['vsv'] , usemask=False)
                    elif 'as' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], np.divide(self.data['vsh'] - self.data['vsv'],self.data['vs'],out=np.zeros_like(self.data['vs']), where= self.data['vs'] != 0.)*100. , usemask=False)
                    elif 'PH-PV' in parameters[ii] or 'ap' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], self.data['vph'] - self.data['vpv'] , usemask=False)
                    elif 'ap' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], np.divide(self.data['vph'] - self.data['vpv'],self.data['vp'],out=np.zeros_like(self.data['vp']), where= self.data['vp'] != 0.)*100. , usemask=False)
                    elif '(SH+SV)*0.5' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], (self.data['vsh'] + self.data['vsv'])/2. , usemask=False)
                    elif '(PH+PV)*0.5' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], (self.data['vph'] + self.data['vpv'])/2. , usemask=False)
                    elif 'dETA/ETA' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], self.data['eta'] , usemask=False)
                    elif 'dRHO/RHO' in parameters[ii]:
                        self.data=append_fields(self.data, parameters[ii], self.data['rho'] , usemask=False)                
                    else:
                        raise NotImplementedError('parameter ',parameters[ii],' is not currently implemented in reference1D.get_custom_parameter')
        else:
            raise ValueError('reference1D object is not allocated')

    def evaluate_at_depth(self,depth_in_km,parameter='vs',interpolation='linear'):   
        '''
        Get the values of a parameter at a given depth
        '''
        values=None
        if isinstance(depth_in_km, (list,tuple,np.ndarray)):
            depth_in_km = np.asarray(depth_in_km)
        elif isinstance(depth_in_km, float):
            depth_in_km = np.asarray([depth_in_km])
        elif isinstance(depth_in_km, int):
            depth_in_km = np.asarray([float(depth_in_km)])
        else:
            raise TypeError('depth_in_km must be list or tuple, not %s' % type(depth_in_km))


        if self.data is not None and self.__nlayers__ > 0:
            if parameter in self.data.dtype.names:
                values = self.data[parameter]
                depth_array = (constants.R - self.data['radius'])/1000. # in km
                # Sort to make interpolation possible
                indx = depth_array.argsort()
                values = griddata(depth_array[indx], values[indx], depth_in_km, method=interpolation)
                if len(depth_in_km)==1: values = values.item()
            else:
                raise ValueError('parameter '+parameter+' not defined in array')
        else:
            raise ValueError('reference1D object is not allocated')
        return values

    def to_cards(self,dir='.',fmt='cards',parameters = ['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta']):
        '''
        Writes a model file that is compatible with MINEOS.
        '''
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            ntotlev = self.__nlayers__
            itopic = self.metadata['discontinuities']['itopic']
            itopoc = self.metadata['discontinuities']['itopoc']
            itopmantle = self.metadata['discontinuities']['itopmantle']
            itopcrust = self.metadata['discontinuities']['itopcrust']
            
            f = open(dir+'/'+model_name+'.'+fmt,'w')
            f.write(model_name+'\n')
            f.write('1 1. 1 1\n')
            line = ff.FortranRecordWriter('(5I5)')
            f.write(line.write([ntotlev,itopic,itopoc,itopmantle,itopcrust])+u'\n')
            line = ff.FortranRecordWriter('(f8.0,3f9.2,2f9.1,2f9.2,f9.5)')

            write = self.data[parameters]
            for i in range(0,len(write)):
                f.write(line.write(write[i])+u'\n')
            f.close()
        else:
            raise ValueError('reference1D object is not allocated')


    def to_TauPmodel(self,dir='.',fmt='tvel'):
        '''
        Writes a model file that is compatible with TauP.
        file format options 'tvel' and 'nd'.

        Note: TauP can't handle zero shear velocity in the ocean layer...
          To work around this, zero values an ocean layer will be written 
          as 1e-4.
        '''
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            f = open(dir+'/'+model_name+'.'+fmt,'w')
            f.write('{} - P\n'.format(model_name))
            f.write('{} - S\n'.format(model_name))

            for i in range(0,len(self.data)):
                f.write('{:2.4f}   {:2.4f}   {:2.4f}    {:2.4f}\n'.format(
                   (self.radius_max - self.data['radius'][::-1][i]) / 1000.0,
                   self.data['vp'][::-1][i] / 1000.0,
                   self.data['vs'][::-1][i] / 1000.0,
                   self.data['rho'][::-1][i] / 1000.0))       
            f.close()
        else:
            raise ValueError('reference1D object is not allocated')

        if self.data['vp'][-1] == 0 or self.data['vs'][-1] == 0:
            raise Warning('zero velocity layer detected at surface ...\n \
                      TauP raytracing may not work')

    def to_axisem(self,dir='.',anelastic=True,anisotropic=True):
        '''
         Write 1D model to be used as an external model in axisem
        '''
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            f = open(dir+'/'+model_name+'.bm','w')
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
        else:
            raise ValueError('reference1D object is not allocated')
            
    def plot(self,figuresize=[7,12],height_ratios=[2, 2, 1],ifshow=True,format='.eps',isotropic=False,zoomdepth=[0.,1000.]):
        """ 
        Plot the cards array in a PREM like plot
        """
        depthkmarr = (constants.R - self.data['radius'])/1000. # in km
        #Set default fontsize for plots
        plots.updatefont(10)
        fig = plt.figure(1, figsize=(figuresize[0],figuresize[1]))
        gs = gridspec.GridSpec(3, 1, height_ratios=height_ratios) 
        fig.patch.set_facecolor('white')
        ax01=plt.subplot(gs[0])
        ax01.plot(depthkmarr,self.data['rho']/1000.,'k')
        ax01.plot(depthkmarr,self.data['vsv']/1000.,'b')
        ax01.plot(depthkmarr,self.data['vsh']/1000.,'b:')
        ax01.plot(depthkmarr,self.data['vpv']/1000.,'r')
        ax01.plot(depthkmarr,self.data['vph']/1000.,'r:')
        mantle=np.where( depthkmarr < 2891.)
        ax01.plot(depthkmarr[mantle],self.data['eta'][mantle],'g')
        ax01.set_xlim([0., constants.R/1000.])        
        ax01.set_ylim([0, 14])
    
        majorLocator = MultipleLocator(2)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(1)
        ax01.yaxis.set_major_locator(majorLocator)
        ax01.yaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax01.yaxis.set_minor_locator(minorLocator)
        
        majorLocator = MultipleLocator(2000)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(1000)
        ax01.xaxis.set_major_locator(majorLocator)
        ax01.xaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax01.xaxis.set_minor_locator(minorLocator)
        ax01.set_ylabel('Velocity (km/sec), density (g/cm'+'$^3$'+') or '+'$\eta$')
        
        for para,color,xloc,yloc in [("$\eta$",'g',1500.,2.),("$V_S$",'b',1500.,7.8),("$V_P$",'r',1500.,13.5),("$\\rho$",'k',1500.,4.5),("$V_P$",'r',4000.,9.2),("$\\rho$",'k',4000.,12.5),("$V_S$",'b',5500.,4.5)]:
            ax01.annotate(para,color=color,
            xy=(3, 1), xycoords='data',
            xytext=(xloc/constants.R/1000., yloc/14.), textcoords='axes fraction',
            horizontalalignment='left', verticalalignment='top')


        ax11=plt.subplot(gs[1])
        depthselect=np.intersect1d(np.where( depthkmarr >= zoomdepth[0]),np.where( depthkmarr <= zoomdepth[1]))
        ax11.plot(depthkmarr[depthselect],self.data['rho'][depthselect]/1000.,'k')
        if isotropic:
            ax11.plot(depthkmarr[depthselect],self.data['vs'][depthselect]/1000.,'b')
        else:
            ax11.plot(depthkmarr[depthselect],self.data['vsv'][depthselect]/1000.,'b')
            ax11.plot(depthkmarr[depthselect],self.data['vsh'][depthselect]/1000.,'b:')
        ax12 = ax11.twinx()
        if isotropic:
            ax12.plot(depthkmarr[depthselect],self.data['vp'][depthselect]/1000.,'r')
        else:
            ax12.plot(depthkmarr[depthselect],self.data['vpv'][depthselect]/1000.,'r')
            ax12.plot(depthkmarr[depthselect],self.data['vph'][depthselect]/1000.,'r:')
        
        ax11.plot(depthkmarr[depthselect],self.data['eta'][depthselect],'g')
        ax11.set_xlim(zoomdepth)
        ax11.set_ylim([0, 7])
        ax12.set_xlim(zoomdepth)
        ax12.set_ylim([-2, 12])        
        ax11.set_ylabel('Shear velocity (km/sec), density (g/cm'+'$^3$'+') or '+'$\eta$')
        ax12.set_ylabel('Compressional velocity (km/sec)')
        for para,color,xloc,yloc in [("$\eta$",'g',150.,1.),("$V_{S}$",'b',150.,4.3),("$V_{P}$",'r',120.,5.5),("$\\rho$",'k',150.,3.8)]:
            ax11.annotate(para,color=color,
            xy=(3, 1), xycoords='data',
            xytext=(xloc/1000., yloc/7.), textcoords='axes fraction',
            horizontalalignment='left', verticalalignment='top')
        ax12.set_yticks(np.arange(6, 14, step=2))
        majorLocator = MultipleLocator(200)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(100)
        ax11.xaxis.set_major_locator(majorLocator)
        ax11.xaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax11.xaxis.set_minor_locator(minorLocator)

        
        ax21=plt.subplot(gs[2], sharex=ax11)
        with np.errstate(divide='ignore', invalid='ignore'): # Ignore warning about dividing by zero
            anisoVs=(self.data['vsh']-self.data['vsv'])*200./(self.data['vsh']+self.data['vsv'])
        anisoVp=(self.data['vph']-self.data['vpv'])*200./(self.data['vph']+self.data['vpv'])
        ax21.plot(depthkmarr[depthselect],anisoVs[depthselect],'b')
        ax21.plot(depthkmarr[depthselect],anisoVp[depthselect],'r')
        ax21.set_ylim([0, 5])        
        ax21.set_xlim(zoomdepth)
        majorLocator = MultipleLocator(1)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(0.5)
        ax21.yaxis.set_major_locator(majorLocator)
        ax21.yaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax21.yaxis.set_minor_locator(minorLocator)
        for para,color,xloc,yloc in [('Q'+'$_{\mu}$','k',400.,2.5),("$a_{S}$",'b',150.,3.7),("$a_{P}$",'r',100.,1.8)]:
            ax21.annotate(para,color=color,
            xy=(3, 1), xycoords='data',
            xytext=(xloc/1000., yloc/4.), textcoords='axes fraction',
            horizontalalignment='left', verticalalignment='top')


        ax22 = ax21.twinx()
        ax22.plot(depthkmarr[depthselect],self.data['Qmu'][depthselect],'k')
        ax21.set_xlabel('Depth (km)')
        ax21.set_ylabel("$V_P$"+' or '+"$V_S$"+' anisotropy (%)')
        ax22.set_ylabel('Shear attenuation Q'+'$_{\mu}$')
        ax22.set_ylim([0, 400])        
        ax22.set_xlim(zoomdepth)
        majorLocator = MultipleLocator(100)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(50)
        ax22.yaxis.set_major_locator(majorLocator)
        ax22.yaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax22.yaxis.set_minor_locator(minorLocator)
        if ifshow: 
            plt.show()
        else:
            plt.savefig(self.name+format)

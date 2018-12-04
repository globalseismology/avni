#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import argparse #parsing arguments
import glob # pattern matching for finding files
import numpy as np #for numerical analysis
from datetime import date  #to give a timestamp to output and compare times
import pdb    #for the debugger pdb.set_trace()
from math import pi
import fortranformat as ff #reading/writing fortran formatted text
from future.utils import native_str
from six import string_types # to check if variable is string using isinstance
from numpy.lib.recfunctions import append_fields # for appending named columns to named numpy arrays
from scipy.interpolate import griddata
import ntpath #Using os.path.split or os.path.basename as others suggest won't work in all cases
import matplotlib.pyplot as plt
from matplotlib import gridspec # Relative size of subplots
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import xarray as xr

####################### IMPORT REM3D LIBRARIES  #######################################
#from . import plots 
#######################################################################################
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
    except IOError:
        raise IOError("File (",filename,") does not exist in the current directory - ",currentdir)

    return epixarr
#####################
# 1D model class

class reference1D(object):
    '''
    A class for 1D reference Earth models used in tomography
    '''

    def __init__(self):
        self.__nlayers__ = None
        self.data = None
        self.metadata = {}
        self.name = None
        self.radius_max = None
    
    def __str__(self):
        if self.data is not None and self.__nlayers__ > 0:
            output = "%s is a one-dimensional model with %s layers and radius up to %s km" % (self.name, self.__nlayers__,self.radius_max/1000.)
        else:
            output = "No model has been read into this reference1D instance yet"
        return output
            
    def read(self,file,fmt='card'):
        '''
        Read a card deck file used in OBANI. Other formats not ready yet
        '''
        if fmt=='card':
            names=['radius','rho','vpv','vsv','Qkappa','Qmu','vph','vsh','eta']
            formats=[np.float for ii in range(len(names))]
            modelarr = np.genfromtxt(file,dtype=None,comments='#',skip_header=3,
            names=names)
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
            self.data=append_fields(self.data, 'A', self.data['rho']*self.data['vph']**2 , usemask=False)
            self.data=append_fields(self.data, 'C', self.data['rho']*self.data['vpv']**2 , usemask=False)
            self.data=append_fields(self.data, 'N', self.data['rho']*self.data['vsh']**2 , usemask=False)
            self.data=append_fields(self.data, 'L', self.data['rho']*self.data['vsv']**2 , usemask=False)
            self.data=append_fields(self.data, 'F', self.data['eta']*(self.data['A']-2.*self.data['L']) , usemask=False)
            self.data=append_fields(self.data, 'kappa', (4.0*(self.data['A']+self.data['F']-self.data['N'])+self.data['C'])/9. , usemask=False)
            self.data=append_fields(self.data, 'mu', (self.data['A']+self.data['C']-2.*self.data['F']+5.*self.data['N']+6.*self.data['L'])/15. , usemask=False)
            self.data=append_fields(self.data, 'vp', np.sqrt(np.divide((self.data['kappa']+4.*self.data['mu']/3.),self.data['rho'])) , usemask=False)
            self.data=append_fields(self.data, 'vs', np.sqrt(np.divide(self.data['mu'],self.data['rho'])) , usemask=False)
            with np.errstate(divide='ignore', invalid='ignore'): # Ignore warning about dividing by zero
                xi = np.power(np.divide(self.data['vsh'],self.data['vsv']),2)
            self.data=append_fields(self.data, 'xi', xi , usemask=False)
            self.data=append_fields(self.data, 'phi', np.power(np.divide(self.data['vpv'],self.data['vph']),2) , usemask=False)
        else:
            raise ValueError('reference1D object is not allocated')

    def get_custom_parameter(self,parameters):
        '''
        Get the arrays of custom parameters defined in various Earth models
        '''
        if self.data is not None and self.__nlayers__ > 0:
            # convert to array for ease of looping
            if isinstance(parameters,string_types): parameters = np.array(parameters) 
            
            for parameter in parameters:
                if 'SH-SV' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['vsh'] - self.data['vsv'] , usemask=False)
                elif 'PH-PV' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['vph'] - self.data['vpv'] , usemask=False)
                elif '(SH+SV)*0.5' in parameter:
                    self.data=append_fields(self.data, parameter, (self.data['vsh'] + self.data['vsv'])/2. , usemask=False)
                elif '(PH+PV)*0.5' in parameter:
                    self.data=append_fields(self.data, parameter, (self.data['vph'] + self.data['vpv'])/2. , usemask=False)
                elif 'dETA/ETA' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['eta'] , usemask=False)
                elif 'dRHO/RHO' in parameter:
                    self.data=append_fields(self.data, parameter, self.data['rho'] , usemask=False)                
                else:
                    raise NotImplementedError('parameter ',parameter,' is not currently implemented in reference1D.get_custom_parameter')
        else:
            raise ValueError('reference1D object is not allocated')

    def evaluate_at_depth(self,depth_in_km,parameter='vs',interpolation='linear'):   
        '''
        Get the values of a parameter at a given depth
        '''
        values=None
        if self.data is not None and self.__nlayers__ > 0:
            if parameter in self.data.dtype.names:
                values = self.data[parameter]
                depth_array = (6371000. - self.data['radius'])/1000. # in km
                # Sort to make interpolation possible
                indx = depth_array.argsort()
                # convert to array for ease of looping
                if isinstance(depth_in_km,float) or isinstance(depth_in_km,int): depth_in_km = np.array(depth_in_km) 
                values = griddata(depth_array[indx], values[indx], depth_in_km, method=interpolation)
            else:
                raise ValueError('parameter '+parameter+' not defined in array')
        else:
            raise ValueError('reference1D object is not allocated')
        return values

    def to_TauPmodel(self,fmt='tvel'):
        '''
        Writes a model file that is compatible with TauP.
        file format options 'tvel' and 'nd'.

        Note: TauP can't handle zero shear velocity in the ocean layer...
        '''
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            f = open(model_name+'.tvel','w')
            f.write(u'{} - P\n'.format(model_name))
            f.write(u'{} - S\n'.format(model_name))

            for i in range(0,len(self.data)):
                f.write(u'{:<10.4f} {:<10.4f} {:<10.4f} {:<10.4f}\n'.format(
                   (self.radius_max - self.data['radius'][::-1][i]) / 1000.0,
                   self.data['vpv'][::-1][i] / 1000.0,
                   self.data['vsv'][::-1][i] / 1000.0,
                   self.data['rho'][::-1][i] / 1000.0))       
            f.close()
        else:
            raise ValueError('reference1D object is not allocated')

        if self.data['vpv'][-1] == 0 or self.data['vsv'][-1] == 0:
            raise Warning('zero velocity layer detected at surface ...\n \
                      TauP raytracing may not work')

    def to_axisem(self,anelastic=True,anisotropic=True):
        '''
         Write 1D model to be used as an external model in axisem
        '''
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            #f = open(model_name+'.bm','w',encoding='utf-8')
            f = open(model_name+'.bm','w')
            n_discon = 0

            if anelastic:
                #f.write(unicode('ANELASTIC     T\n'))
                f.write(u'ANELASTIC     T\n')
            else:
                #f.write(unicode('ANELASTIC     F\n'))
                f.write(u'ANELASTIC     F\n')

            if anisotropic:
                #f.write(unicode('ANISOTROPIC     T\n'))
                f.write(u'ANISOTROPIC     T\n')
            else:
                #f.write(unicode('ANISOTROPIC     F\n'))
                f.write(u'ANISOTROPIC     F\n')

            #f.write(unicode('UNITS      m\n'))
            f.write(u'UNITS      m\n')

            if anisotropic:
                #f.write(unicode('COLUMNS   radius    rho    vpv    vsv    qka    qmu    vph    vsh    eta\n'))
                f.write(u'COLUMNS   radius    rho    vpv    vsv    qka    qmu    vph    vsh    eta\n')

            for i in range(0,len(self.data)):
                f.write(u'{:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f}\n'.format(
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
                    f.write(u'#    Discontinuity {}, depth {:6.2f} km\n'.format(n_discon,depth_here))
        else:
            raise ValueError('reference1D object is not allocated')
            
    def plot(self,figuresize=[7,12],height_ratios=[2, 2, 1],ifshow=True,outfile='.png'):
        """ 
        Plot the cards array in a PREM like plot
        """
        depthkmarr = (6371000. - self.data['radius'])/1000. # in km
        #Set default fontsize for plots
        plots.updatefont(None,10)
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
        ax01.set_xlim([0., 6371.])
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
            xytext=(xloc/6371., yloc/14.), textcoords='axes fraction',
            horizontalalignment='left', verticalalignment='top')


        ax11=plt.subplot(gs[1])
        top1000km=np.where( depthkmarr < 1000.)
        ax11.plot(depthkmarr[top1000km],self.data['rho'][top1000km]/1000.,'k')
        ax11.plot(depthkmarr[top1000km],self.data['vsv'][top1000km]/1000.,'b')
        ax11.plot(depthkmarr[top1000km],self.data['vsh'][top1000km]/1000.,'b:')
        ax12 = ax11.twinx()
        ax12.plot(depthkmarr[top1000km],self.data['vpv'][top1000km]/1000.,'r')
        ax12.plot(depthkmarr[top1000km],self.data['vph'][top1000km]/1000.,'r:')
        ax11.plot(depthkmarr[top1000km],self.data['eta'][top1000km],'g')
        ax11.set_xlim([0., 1000.])
        ax11.set_ylim([0, 7])
        ax12.set_xlim([0., 1000.])
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
        ax21.plot(depthkmarr[top1000km],anisoVs[top1000km],'b')
        ax21.plot(depthkmarr[top1000km],anisoVp[top1000km],'r')
        ax21.set_ylim([0, 4])        
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
        ax22.plot(depthkmarr[top1000km],self.data['Qmu'][top1000km],'k')
        ax21.set_xlabel('Depth (km)')
        ax21.set_ylabel("$V_P$"+' or '+"$V_S$"+' anisotropy (%)')
        ax22.set_ylabel('Shear attenuation Q'+'$_{\mu}$')
        ax22.set_ylim([0, 400])        
        majorLocator = MultipleLocator(100)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(50)
        ax22.yaxis.set_major_locator(majorLocator)
        ax22.yaxis.set_major_formatter(majorFormatter)
        # for the minor ticks, use no labels; default NullFormatter
        ax22.yaxis.set_minor_locator(minorLocator)
        if ifshow: plt.show()
        plt.savefig(self.name+outfile)

def epix2ascii(epix_dir,output_file,ref_model,mod_desc,n_hpar=1):
    '''
    write a rem3d formatted ascii file from a directory containing epix files 

    Input parameters:
    ----------------
  
    epix_dir: path to directory containing epix layer files
    
    output_file: name of rem3d format output file
    
    ref_model: the 1D reference model used for the specified tomography model
    
    mod_desc: tomographic model parameter (e.g. "(SH+SV)*0.5")
    
    n_hpar:number of horizontal parameterizations (currently only handles 
           models with 1 horizontal parameterization, and with a constant pixel width)
           
    '''
    f_out = open(output_file,'w')
    cwd = os.getcwd()

    #find layer info
    globpath=cwd+'/'+epix_dir+'/*'
    filenames=glob.glob(globpath) 
    layer_start = []
    layer_end = []

    for filename in filenames:
        print(filename)
        start_ = filename.split('.epix')[-2].split('_')[-2]
        end_ = filename.split('.epix')[-2].split('_')[-1]
        layer_start.append(float(start_))
        layer_end.append(float(end_))

    layer_start.sort()
    layer_end.sort()
    zmin = np.min(layer_start)
    zmax = np.max(layer_end)
    dz = layer_end[1] - layer_start[1]
    depths = np.arange(dz,2900.0+dz,dz)
    depths[0] = zmin
    depths[-1] = zmax

    #write header
    f_out.write(u'REFERENCE MODEL: {} \n'.format(ref_model))
    f_out.write(u'KERNEL SET: {}\n'.format(kernel_set))

    #write radial kernels
    n_rad_kernels = len(depths)-1
    f_out.write(u'RADIAL STRUCTURE KERNELS: {}\n'.format(n_rad_kernels))
    for i in range(0,n_rad_kernels):
        f_out.write(u'DESC  {:3.0f}: {}, {:1.1f} - {:1.1f} km\n'.format(i+1,
        mod_desc,depths[i],depths[i+1]))

    #write horizontal parameterization(s)
    epixarr = np.loadtxt(filenames[0])
    lats = epixarr[:,0]
    lons = epixarr[:,1]
    px_w = epixarr[:,2] 
    if len(np.unique(pixel_width)) is not 1: 
        raise ValueError('epix2ascii can only handle constant pixel width')
    pixel_width = px_w[0]
    kernel_set = 'BOX{:1.0f}+I1D'.format(dz)
    #kernel_set = 'BOX{:1.0f}_PIX{:1.1f}'.format(dz,pixel_width)

    f_out.write(u'HORIZONTAL PARAMETERIZATIONS: {}\n'.format(n_hpar))
    for i in range(0,n_hpar):
        f_out.write(u'HPAR   {}: PIXELS,  {:1.1f} x {:1.1f}\n'.format(                               i+1,pixel_width,pixel_width))
        shape = (int(180.0/pixel_width),int(360.0/pixel_width)) 
        lats_arr = np.reshape(lats,shape,order='F')
        lons_arr = np.reshape(lons,shape,order='F')
        px_w_arr = np.reshape(px_w,shape,order='F')

        lats = lats_arr.flatten()
        lons = lons_arr.flatten()
        px_w = px_w_arr.flatten()

    for j in range(0,len(lons)):
        if lons[i] > 180.0:
            lons[i] -= 360.0
            f_out.write(u'{:5.1f} {:5.1f} {:5.1f}\n'.format(lons[i],
                    lats[i], px_w[i]))

    #write model coefficients
    line = ff.FortranRecordWriter('(6E12.4)')
    for i in range(0,len(depths)-1):
        print('writing coefficients for layer ', i+1)
        epix_name='*{:1.1f}_{:1.1f}.epix'.format(depths[i],depths[i+1])
        epix_glob = glob.glob(cwd+'/'+epix_dir+'/'+epix_name)
        print(epix_glob)
        epixarr = np.loadtxt(epix_glob[0])
        coefs = epixarr[:,3]
        coefs_arr = np.reshape(coefs,shape,order='F')
        coefs = coefs_arr.flatten()
        f_out.write(u'STRU  {:3.0f}:  {:1.0f}\n'.format(i+1,pixel_width))
        f_out.write(line.write(coefs)+u'\n')

def ascii2xarray(ascii_file,save_netcdf=False,outfile=None):
    '''
    write a rem3d formatted ascii file from a directory containing epix files 

    Input parameters:
    ----------------
  
    ascii_file: path to rem3d format output file
    
    save_netcdf: save netcdf format
    
    outfile: output netcdf file
    
    '''

    with open(ascii_file) as f:
    #read header
        for i, line in enumerate(f):
            if i == 0:
                ref_model = line.split('REFERENCE MODEL:')[1].strip()
            elif i == 1:
                krnl_set = line.split('KERNEL SET:')[1].strip()
            elif i == 2:
                nrad_krnl = line.split('RADIAL STRUCTURE KERNELS:')[1].strip()
                nrad_krnl = int(nrad_krnl)
            else:
                break

        #read radial kernels
        rkrnl_start = np.zeros(nrad_krnl)
        rkrnl_end = np.zeros(nrad_krnl)
        #TODO implement other radial kernel options
        if krnl_set[0:3] == 'BOX': krnl_wdth = float(krnl_set[3:].split('+')[0])
        npts_dep = nrad_krnl

        for i, line in enumerate(f):
            if i <= nrad_krnl-1:
                mod_par = line.strip().split()[2].split(',')[0]
                rkrnl_start[i] = float(line.strip().split()[3])
                rkrnl_end[i] = float(line.strip().split()[5])
            if i == nrad_krnl-1: break

        #read horizontal parameterization info
        hpar_list = []
        #TODO implement parameterizations besides pixels
        for i, line in enumerate(f):
            if i == 0:
                nhpar = int(line.strip().split()[-1])
                print('NHPAR', nhpar)
            elif i== 1:
                hpar_type = line.strip().split()[2].split(',')[0]
                hpar_name = line.split(':')[1].strip()
                hpar_list.append(hpar_name)

                if hpar_type == 'PIXELS':
                    lons = []
                    lats = []
                    pxwd = []
                    pxw_lon = float(line.strip().split()[3])
                    pxw_lat = float(line.strip().split()[5])
                break

        for j in range(0,nhpar):
            for i, line in enumerate(f):

                if line.strip().split()[0] == 'HPAR':
                    par_type = line.strip()[2].split(',')[0]
                    #if par_type == 'PIXELS':
                    #    pxw_lon = float(line.strip()[3])
                    #    pxw_lat = float(line.strip()[4])
                    #    npts_lon = int(360. / pxw_lon)
                    #        npts_lat = int(180. / pxw_lat)
                    break

                elif line.strip().split()[0] == 'STRU':
                    i_layer = int(line.strip().split()[1].split(':')[0])
                    i_hpar = int(line.strip().split()[2])
                    break

                else:
                    lons.append(float(line.strip().split()[0]))
                    lats.append(float(line.strip().split()[1]))
                    pxwd.append(float(line.strip().split()[2]))

        #read structure coefficients block
        layer_coefs = []
        layer_dict = {}
        for i, line in enumerate(f):

            if line.strip().split()[0] == 'STRU':
                layer_dict[i_layer] = layer_coefs
                i_layer = int(line.strip().split()[1].split(':')[0])
                i_hpar = int(line.strip().split()[2])
                layer_coefs = []
            else:
                for coef in line.strip().split(): 
                    layer_coefs.append(float(coef))

        layer_dict[i_layer] = layer_coefs #add layer layer to dict
    
    #create dims arrays
    lon = np.arange((pxw_lon/2.),360.,pxw_lon)
    for i in range(0,len(lon)):
        if lon[i] >= 180.0:
            lon[i] -= 360.0
    lat = np.arange(-90.+(pxw_lat/2.),90,pxw_lat)
    dep = (rkrnl_start + rkrnl_end) / 2.

    #create xarray 
    data_array = xr.DataArray(np.zeros((len(dep),len(lon),len(lat))),
                              dims = ['depth','lon','lat'])

    data_array['lon'] = lon
    data_array['lat'] = lat
    data_array['depth'] = dep

    #fill in data array
    for i,layer in enumerate(layer_dict):
        data_array[i,:,:] = np.reshape(layer_dict[layer],
                                (len(lon),len(lat)),order='F')
    print('done reading model')
        
    if save_netcdf: data_array.to_netcdf(outfile)

    return data_array

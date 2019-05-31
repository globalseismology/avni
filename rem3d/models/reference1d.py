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
import re
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
        self.__nlayers__ = None
        self.data = None
        self.metadata = {}
        # assume that information about the native parameterization is not available
        # this is typical for a card deck file
        for field in ['model','ref_period','parameters']: self.metadata[field] = None
        self.name = None
        self.radius_max = None
        if file is not None:
            self.read(file)
            if self.data is not None:
                self.get_Love_elastic()
                self.get_discontinuity()
                self.get_mineralogical()

    def __str__(self):
        if self.data is not None and self.__nlayers__ > 0:
            output = "%s is a one-dimensional model with %s layers and radius up to %s km" % (self.name, self.__nlayers__,self.radius_max/1000.)
        else:
            output = "No model has been read into this reference1D instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self.name})'.format(self=self)

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
            self.read_mineos_cards(file)
        except:
            var1 = traceback.format_exc()
            try:
                self.read_bases_coefficients(file)
            except:
                var2 = traceback.format_exc()
                print('############    Tried reading as mineos cards   ############')
                print(var1)
                print('############    Tried reading as bases coefficients   ############')
                print(traceback.format_exc())
                raise NotImplementedError('model format is not currently implemented in reference1D.read')

    def read_bases_coefficients(self,file):
        # Use tools.eval_polynomial and tools.eval_splrem
        coef_names=['bottom','top','spline','constant','linear','quadratic','cubic']
        # function used to parse line with key word from rx_dict
        def _parse_line(line,rx_dict):
            """
            Do a regex search against all defined regexes and
            return the key and match result of the first matching regex
            """
            for key, rx in rx_dict.items():
                match = rx.search(line)
                if match:
                    return key, match
            # if there are no matches
            return None, None

        # regex dict for common metadata info
        rx_dict_common = {
            'model': re.compile(r'EARTH MODEL\s*:\s*(?P<model>.*)\n'),
            'ref_period': re.compile(r'REFERENCE PERIOD\s*:\s*(?P<ref_period>.*)\n'),
            'norm_radius': re.compile(r'NORMALIZING RADIUS\s*:\s*(?P<norm_radius>.*)\n'),
            'parameters': re.compile(r'#PARAMETERS\s*:\s*(?P<parameters>.*)\n'),
            'num_region': re.compile(r'NUMBER OF REGIONS\s*:\s*(?P<num_region>.*)\n'),
            'regions': re.compile(r'REGION\s*:\s*(?P<regions>.*)\n'),
        }

        # Check if it is the first file read for the model
        if self.metadata['model'] == None:
            # make parameterization 2D lists of dicts,first dimension associated with file
            # number, here we assume the parameterization within each file are the same
            self.metadata['parameterization'] = [[]]
            # make parameters list of dicts, PARAMETERS should not be the same
            self.metadata['parameters'] = {}
            regions = []
            # loop through lines in file, extract info
            with open(file,'r') as f:
                line = f.readline()
                while line:
                # at each line check for a match with a regex
                    key, match = _parse_line(line,rx_dict_common)
                    if key == 'model':
                        self.metadata['model'] = match.group('model')
                        self.name = self.metadata['model']
                    if key == 'ref_period':
                        ref_temp = match.group('ref_period')
                        self.metadata['ref_period'] = float(ref_temp)
                    if key == 'norm_radius':
                        rad_temp = match.group('norm_radius')
                        self.metadata['norm_radius'] = float(rad_temp)
                    if key == 'parameters':
                        para_list = match.group('parameters').split()
                        self.metadata['parameter_list']=para_list
                    if key == 'num_region':
                        nr_temp = match.group('num_region')
                        num_region = int(nr_temp)
                    if key == 'regions':
                        regions.append(match.group('regions').strip())
                    line = f.readline()
        else:
            # append a new parameterization
            self.metadata['parameterization'].append([])
            regions = []
            with open(file,'r') as f:
                line = f.readline()
                while line:
                # at each line check for a match with a regex
                    key, match = _parse_line(line,rx_dict_common)
                    # Check if model name is the same
                    if key == 'model':
                        if self.name != match.group('model'):
                            raise ValueError('model names should match between input files')
                    if key == 'parameters':
                        para_list = match.group('parameters').split()
                        self.metadata['parameter_list'].append(para_list)
                    if key == 'num_region':
                        nr_temp = match.group('num_region')
                        num_region = int(nr_temp)
                    if key == 'regions':
                        regions.append(match.group('regions'))
                    line = f.readline()
        # Now start read parameterization from the file
        rx_dict_para = {
            'bot_dep': re.compile(r'BOT DEPTH\s*:\s*(?P<bot_dep>.*)\n'),
            'top_dep': re.compile(r'TOP DEPTH\s*:\s*(?P<top_dep>.*)\n'),
            'bot_rad': re.compile(r'BOT RADIUS\s*:\s*(?P<bot_rad>.*)\n'),
            'top_rad': re.compile(r'TOP RADIUS*:\s*(?P<top_rad>.*)\n'),
            'lvl': re.compile(r'LEVELS\s*:\s*(?P<lvl>.*)\n'),
            'poly': re.compile(r'POLYNOMIAL\s*:\s*(?P<poly>.*)\n'),
        }
        bot_deps = []
        top_deps = []
        bot_rads = []
        top_rads = []
        levels = []
        polys = []
        #radius_flag = 0 #The depth info could be saved as radius or depth, use this to mark
        with open(file,'r') as f:
            line = f.readline()
            while line:
            # at each line check for a match with a regex
                key, match = _parse_line(line,rx_dict_para)
                if key == 'poly':
                    polys.append(match.group('poly'))
                if key == 'lvl':
                    levels.append(int(match.group('lvl')))
                if key == 'bot_dep':
                    bd_temp=match.group('bot_dep')
                    bot_rads.append(self.metadata['normalizing_radius']-float(bd_temp))
                if key == 'top_dep':
                    td_temp=match.group('top_dep')
                    top_rads.append(self.metadata['normalizing_radius']-float(td_temp))
                if key == 'bot_rad':
                    br_temp=match.group('bot_rad')
                    bot_rads.append(float(br_temp))
                    radius_flag = 1
                if key == 'top_rad':
                    tr_temp=match.group('top_rad')
                    top_rads.append(float(tr_temp))
                line = f.readline()
        if radius_flag == 0:
            bot_rads = self.metadata['norm_radius']-np.array(bot_deps)
            top_rads = self.metadata['norm_radius']-np.array(top_deps)
        elif radius_flag == 1:
            bot_rads = np.array(bot_rads)
            top_rads = np.array(top_rads)
        # assign the num of regions to the n th parameterization
        iparam = len(self.metadata['parameterization']) - 1
        self.metadata['parameterization'][iparam] = {'num_regions':num_region,'filename':file,'description':'read from '+file}
        for idx,(region,poly,level,bot_rad,top_rad) in enumerate(zip(regions,polys,levels,bot_rads,top_rads)):
            self.metadata['parameterization'][iparam].update({region:{'polynomial':poly,'levels':level,'top_radius':top_rad,'bottom_radius':bot_rad}})

        # Now start to read parameter coefficient from the file
        # regex for model coefficient info
        rx_dict_coef = {
            'regions': re.compile(r'REGION\s*:\s*(?P<regions>.*)\n'),
            'poly': re.compile(r'POLYNOMIAL\s*:\s*(?P<poly>.*)\n'),
            'comment':re.compile(r'#--*\n')
        }
        for para in para_list:
            self.metadata['parameters'].update({para:{'param_index':iparam}})
        with open(file,'r') as f:
            line = f.readline()
            while line:
            # at each line check for a match with a regex
                key, match = _parse_line(line,rx_dict_coef)
                if key == 'regions':
                    reg_temp = (match.group('regions').strip())
                    for para in para_list: self.metadata['parameters'][para].update({reg_temp:{}})
                if key == 'poly':
                    line=f.readline()
                    while line:
                        key, match = _parse_line(line,rx_dict_coef)
                        if key == 'comment': break
                        att,coef = line.split(":",1)
                        coef = np.array(coef.split())
                        for idx, para in enumerate(para_list):
                            self.metadata['parameters'][para][reg_temp].update({att.strip():coef[idx].astype(float)})
                        line = f.readline()
                line = f.readline()

    def coefficients_to_cards(self):
        """evaluates bases coefficients at prescribed depth levels to get a card deck file
        """
        pdb.set_trace()

    def read_mineos_cards(self,file):
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
        self.__nlayers__ = len(modelarr['radius'])
        # Create data array
        PA_ = pint.PintArray
        modelarr_['depth'] = PA_((constants.R.magnitude - modelarr_['radius'].pint.to(constants.R.units).data).tolist(), dtype = constants.R.units)
        self.data = modelarr_
        self.radius_max = max(self.data['radius']).magnitude

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
        if self.data is not None and self.__nlayers__ > 0:
            if constants.planetpreferred == 'Earth':
                file = self.metadata['filename']
                layers = self.__nlayers__
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

        for field in ['delta','average']: disc[field] = self.data.copy().drop(range(len(np.unique(disc_depths)),len(self.data)))

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
        if self.data is not None and self.__nlayers__ > 0:
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
        radius_in_km = self.metadata['normalizing_radius'] - depth_in_km
        # detailed information about the native parameterization which went into the
        # inversion is available
        if self.metadata['parameters'] is not None:
        # CHAO FILL THIS UP BY USING 
            param_indx = self.metadata['parameters'][parameter.upper()]['param_index']
            for region in self.metadata['parameterization'][param_indx]:
                if region['top_radius']  is not None:
                    dep_flag = (region['top_radius'] - radius_in_km)*(region['bottom_radius']-radius_in_km)
                    if dep_flag <= 0:
                        target_region = region
                        break
            if target_region is None:
                #raise error
            else:
                rnorm = self.metadata['normalizing_radius']
                types = self.metadata['parameterization'][param_indx][target_region]['polynomial']
                radius_range = [self.metadata['parameterization'][param_indx][target_region]['bottom_radius'],
                                self.metadata['parameterization'][param_indx][target_region]['top_radius']]
            # USE this appropriuately
            vercof,dvercof = rem3d.tools.bases.eval_polynomial(radius_in_km,radius_range ,rnorm,types)

        # If detailed metadata regarding the basis parameterization is not available
        # interpolated based on the card deck file
        else:
            if self.data is not None:
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
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            ntotlev = self.__nlayers__
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
        if self.data is not None and self.__nlayers__ > 0:
            model_name = self.name
            f = open(directory+'/'+model_name+'.'+fmt,'w')
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

    def to_axisem(self,directory='.',anelastic=True,anisotropic=True):
        '''
         Write 1D model to be used as an external model in axisem
        '''
        if self.data is not None and self.__nlayers__ > 0:
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
                    depth_here = (self.radius_max - self.data['radius'][::-1][i]) / 1000.0
                    n_discon += 1
                    f.write('#    Discontinuity {}, depth {:6.2f} km\n'.format(n_discon,depth_here))
        else:
            raise ValueError('reference1D object is not allocated')
#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets
in the standard AVNI format."""

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
from scipy.interpolate import griddata,splrep,splev
import copy
from collections import Counter
import traceback
import pandas as pd
import pint
import re
import ntpath
import uuid
import contextlib
import warnings
import os

if sys.version_info[0] >= 3: unicode = str

####################### IMPORT AVNI LIBRARIES  #######################################
from .. import constants
from .. import tools
from avni.f2py import getbullen
#######################################################################################
# 1D model class

class Reference1D(object):
    '''
    A class for 1D reference Earth models used in tomography
    '''
    #########################       magic       ##########################

    def __init__(self,file=None):
        self.data = None
        self.metadata = {}
        # assume that information about the native parameterization is not available
        # this is typical for a card deck file
        for field in ['ref_period','parameters']: self.metadata[field] = None
        self._name = None
        self._radius_max = None
        self._nlayers = None
        if file is not None:
            self.read(file)
            self.derive()

    def __str__(self):
        if self.data is not None and self._nlayers > 0:
            output = "%s is a one-dimensional model with %s layers and radius up to %s km" % (self._name, self._nlayers,self._radius_max/1000.)
        else:
            output = "No model has been read into this reference1D instance yet"
        return output

    def __repr__(self):
        return '{self.__class__.__name__}({self._name})'.format(self=self)

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
            setattr(result, k, copy.deepcopy(v, memo))
        return result

    def __add__(self, other):
        raise NotImplementedError('method to add 1D instances on top of each other')

    #########################       decorators       ##########################

    @property
    def name(self):
        return self._name

    #########################       methods       #############################

    def derive(self):
        if self.data is not None and self._nlayers > 0:
            self.get_Love_elastic()
            self.get_discontinuity()
            self.get_custom_parameter(['as','ap'])

    def read(self,file: str):
        '''
        Read a card deck file used in OBANI. Other formats not ready yet
        '''
        if os.path.isfile(file+'.elas') and os.path.isfile(file+'.anelas'):
            warnings.warn('Reading seperate elas and anelas files as Reference1D. This is preferred since derivatives needed for mineralogical parameters have lower numerical precisions.')
            listfile = [file+'.elas', file+'.anelas']
            for file in listfile: self.read_bases_coefficients(file)
        elif os.path.isfile(file):
            listfile = [file]
        else:
            raise ValueError('Input Reference1D file or seperate constituent files do not exist : '+file)

        try:
            for file in listfile: self.read_bases_coefficients(file)
        except:
            var1 = traceback.format_exc()
            try:
                self.read_mineos_cards(file)
                warnings.warn('Reading as a card file into Reference1D')
            except:
                var2 = traceback.format_exc()
                print('############    Tried reading as basis coefficients   ############')
                print(var1)
                print('############    Tried reading as mineos cards   ############')
                print(var2)
                raise NotImplementedError('model format is not currently implemented in reference1D.read')

    def plot(self):
        from ..plots import plotreference1d
        plotreference1d(self)

    def read_bases_coefficients(self,file: str):
        # regex dict for common metadata info
        rx_dict_common = {
            'model': re.compile(r'EARTH MODEL\s*:\s*(?P<model>.*)\n'),
            'ref_period': re.compile(r'REFERENCE PERIOD\s*:\s*(?P<ref_period>.*)\n'),
            'norm_radius': re.compile(r'NORMALIZING RADIUS\s*:\s*(?P<norm_radius>.*)\n'),
            'parameters': re.compile(r'#PARAMETERS\s*:\s*(?P<parameters>.*)\n'),
            'num_region': re.compile(r'NUMBER OF REGIONS\s*:\s*(?P<num_region>.*)\n'),
            'regions': re.compile(r'REGION\s*:\s*(?P<regions>.*)\n'),
            'reg_codes': re.compile(r'REG CODE\s*:\s*(?P<reg_codes>.*)\n'),
            'units': re.compile(r'#UNITS\s*:\s*(?P<units>.*)\n'),
        }

        # Check if it is the first file read for the model
        if self._name == None:
            # make parameterization 2D lists of dicts,first dimension associated with file
            # number, here we assume the parameterization within each file are the same
            self.metadata['parameterization'] = [[]]
            # make parameters list of dicts, PARAMETERS should not be the same
            self.metadata['attributes'] = {}
            regions = []; reg_codes = []
            # loop through lines in file, extract info
            with open(file,'r') as f:
                line = f.readline()
                while line:
                # at each line check for a match with a regex
                    key, match = tools.parse_line(line,rx_dict_common)
                    if key == 'model':
                        self._name = match.group('model')
                    if key == 'ref_period':
                        ref_temp = match.group('ref_period')
                        self.metadata['ref_period'] = float(ref_temp)
                    if key == 'norm_radius':
                        rad_temp = match.group('norm_radius')
                        self.metadata['norm_radius'] = float(rad_temp)
                        if self.metadata['norm_radius'] != constants.R.to('km').magnitude:
                            raise AssertionError('norm_radius '+str(self.metadata['norm_radius'])+' is consistent with avni constant '+str(constants.R.to('km').magnitude)+'. Reinitialize avni with tools.common.getplanetconstants.')
                    if key == 'parameters':
                        para_list = match.group('parameters').lower().split()
                        self.metadata['parameters']=para_list
                    if key == 'units':
                        unit_list = match.group('units').split()
                        self.metadata['units']=unit_list
                    if key == 'num_region':
                        nr_temp = match.group('num_region')
                        num_region = int(nr_temp)
                    if key == 'regions':
                        regions.append(match.group('regions').strip().lower())
                    if key == 'reg_codes':
                        reg_codes.append(match.group('reg_codes').strip().upper())
                    line = f.readline()
        else:
            # append a new parameterization
            self.metadata['parameterization'].append([])
            regions = []; reg_codes = []
            with open(file,'r') as f:
                line = f.readline()
                while line:
                # at each line check for a match with a regex
                    key, match = tools.parse_line(line,rx_dict_common)
                    # Check if model name is the same
                    if key == 'model':
                        if self._name != match.group('model'):
                            raise ValueError('model names should match between input files')
                    if key == 'ref_period':
                        ref_temp2 = match.group('ref_period')
                        if self.metadata['ref_period'] != float(ref_temp2):
                            print('reference period not consistent!')
                    if key == 'norm_radius':
                        rad_temp2 = match.group('norm_radius')
                        if self.metadata['norm_radius'] != float(rad_temp2):
                            print('normalizing period not consistent!')
                        if self.metadata['norm_radius'] != constants.R.to('km').magnitude:
                            raise AssertionError('norm_radius '+str(self.metadata['norm_radius'])+' is consistent with avni constant '+str(constants.R.to('km').magnitude)+'. Reinitialize avni with tools.common.getplanetconstants.')
                    if key == 'parameters':
                        para_list = match.group('parameters').lower().split()
                        for para in para_list:
                            self.metadata['parameters'].append(para)
                    if key == 'units':
                        unit_list = match.group('units').split()
                        for unit in unit_list:
                            self.metadata['units'].append(unit)
                    if key == 'num_region':
                        nr_temp = match.group('num_region')
                        num_region = int(nr_temp)
                    if key == 'regions':
                        regions.append(match.group('regions').strip().lower())
                    if key == 'reg_codes':
                        reg_codes.append(match.group('reg_codes').strip().upper())
                    line = f.readline()
        # Now start read parameterization from the file
        rx_dict_para = {
            'bot_dep': re.compile(r'BOT DEPTH\s*:\s*(?P<bot_dep>.*)\n'),
            'top_dep': re.compile(r'TOP DEPTH\s*:\s*(?P<top_dep>.*)\n'),
            'bot_rad': re.compile(r'BOT RADIUS\s*:\s*(?P<bot_rad>.*)\n'),
            'top_rad': re.compile(r'TOP RADIUS*:\s*(?P<top_rad>.*)\n'),
            'lvl': re.compile(r'LEVELS\s*:\s*(?P<lvl>.*)\n'),
            'poly': re.compile(r'POLYNOMIAL\s*:\s*(?P<poly>.*)\n'),
            'spln': re.compile(r'SPLINEPNTS\s*:\s*(?P<spln>.*)\n'),
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
                key, match = tools.parse_line(line,rx_dict_para)
                if key == 'poly':
                    polys.append(match.group('poly'))
                if key == 'spln':
                    polys.append('SPLINEPNTS: ' + match.group('spln'))
                if key == 'lvl':
                    levels.append(int(match.group('lvl')))
                if key == 'bot_dep':
                    bd_temp=match.group('bot_dep')
                    bot_rads.append(self.metadata['norm_radius']-float(bd_temp))
                    bot_deps.append(float(bd_temp))
                if key == 'top_dep':
                    td_temp=match.group('top_dep')
                    top_rads.append(self.metadata['norm_radius']-float(td_temp))
                    top_deps.append(float(td_temp))
                if key == 'bot_rad':
                    br_temp=match.group('bot_rad')
                    bot_rads.append(float(br_temp))
                    bot_deps.append(self.metadata['norm_radius']-float(br_temp))
                    radius_flag = 1
                if key == 'top_rad':
                    tr_temp=match.group('top_rad')
                    top_rads.append(float(tr_temp))
                    top_deps.append(self.metadata['norm_radius']-float(tr_temp))
                line = f.readline()
            bot_rads = np.array(bot_rads); top_rads = np.array(top_rads)
            bot_deps = np.array(bot_deps); top_deps = np.array(top_deps)
        # assign the num of regions to the n th parameterization
        iparam = len(self.metadata['parameterization']) - 1
        self.metadata['parameterization'][iparam] = {'num_regions':num_region,'filename':file,'description':'read from '+file}
        for idx,(region,reg_code,poly,level,bot_rad,top_rad,bot_dep,top_dep) in enumerate(zip(regions,reg_codes,polys,levels,bot_rads,top_rads,bot_deps,top_deps)):
            self.metadata['parameterization'][iparam].update({region:{'polynomial':poly,'reg_codes':reg_code,'levels':level,'top_radius':top_rad,'bottom_radius':bot_rad,'top_depth':top_dep,'bottom_depth':bot_dep}})

        # Now start to read parameter coefficient from the file
        # regex for model coefficient info
        rx_dict_coef = {
            'regions': re.compile(r'REGION\s*:\s*(?P<regions>.*)\n'),
            'poly': re.compile(r'(POLYNOMIAL|SPLINEPNTS)\s*:\s*(?P<poly>.*)\n'),
            'comment':re.compile(r'#--*\n')
        }
        for para in para_list:
            self.metadata['attributes'].update({para:{'param_index':iparam}})
        with open(file,'r') as f:
            line = f.readline()
            while line:
            # at each line check for a match with a regex
                key, match = tools.parse_line(line,rx_dict_coef)
                if key == 'regions':
                    reg_temp = (match.group('regions').strip().lower())
                    for para in para_list: self.metadata['attributes'][para].update({reg_temp:{}})
                if key == 'poly':
                    line=f.readline()
                    while line:
                        key, match = tools.parse_line(line,rx_dict_coef)
                        if key == 'comment': break
                        att,coef = line.split(":",1)
                        coef = np.array(coef.split())
                        for idx, para in enumerate(para_list):
                            self.metadata['attributes'][para][reg_temp].update({att.strip():coef[idx].astype(float)})
                        line = f.readline()
                line = f.readline()
        self.metadata['filename'] = file

    def write_bases_coefficients(self,file: str, verbose: bool = False):
        parameterization = self.metadata['parameterization']
        attributes = self.metadata['attributes']
        parameters = self.metadata['parameters']
        units = self.metadata['units']
        npar = len(parameterization)
        if npar > 2 : raise NotImplementedError('More than two parameterizations not availabe for write_bases_coefficients')

        for ipar in range(npar):
            listvar = [key for key in attributes if attributes[key]['param_index']==ipar]
            listunits = [key for ii,key in enumerate(units) if parameters[ii] in listvar]
            if npar > 1:
                if 'vsh' in listvar or 'vs' in listvar:
                    outfile = file+'.elas'
                else:
                    outfile = file+'.anelas'
            else:
                outfile = file

            # Write Header
            f = open(outfile,'w')
            f.write('EARTH MODEL        :{}\n'.format(self._name))
            f.write('REFERENCE PERIOD   :{}\n'.format(self.metadata['ref_period']))
            f.write('NORMALIZING RADIUS :{}\n'.format(self.metadata['norm_radius']))
            f.write('NUMBER OF REGIONS  :{}\n'.format(parameterization[ipar]['num_regions']))

            f.write('#-----------------------------------------------\n')
            strout1 = '#PARAMETERS:'
            strout2 = '#UNITS:     '
            for ii,temp in enumerate(listunits):
                if len(temp) > 10:
                    strout1 += '{:<14}'.format(listvar[ii].strip().upper())
                    strout2 += '{:<14}'.format(temp)
                else:
                    strout1 += '{:<9}'.format(listvar[ii].strip().upper())
                    strout2 += '{:<9}'.format(temp)
            f.write(strout1+'\n')
            f.write(strout2+'\n')
            f.write('#-----------------------------------------------\n')

            # Loop over regions
            listreg = [key for key in parameterization[ipar].keys() if type(parameterization[ipar][key])==dict]
            for reg in listreg:
                f.write('REGION    :{}\n'.format(reg.strip().upper()))
                meta = parameterization[ipar][reg]
                f.write('REG CODE  :{}\n'.format(meta['reg_codes'].strip().upper()))
                if meta['bottom_depth'] > 771:
                    f.write('BOT RADIUS:{}\n'.format(meta['bottom_radius']))
                else:
                    f.write('BOT DEPTH :{}\n'.format(meta['bottom_depth']))
                if meta['top_depth'] > 771:
                    f.write('TOP RADIUS:{}\n'.format(meta['top_radius']))
                else:
                    f.write('TOP DEPTH :{}\n'.format(meta['top_depth']))
                f.write('LEVELS    :{}\n'.format(meta['levels']))
                if meta['polynomial'].lower().startswith('splinepnts'):
                    f.write('{}\n'.format(meta['polynomial']))
                else:
                    f.write('POLYNOMIAL:{}\n'.format(meta['polynomial']))

                # Loop over variables
                basis = attributes[listvar[0]][reg].keys()
                for bas in basis:
                    strout = '{:<10}:'.format(bas.strip().upper())
                    for var in listvar:
                        if var.strip().lower() in ['qkappa','qmu']:
                            temp = '{:7.1f}'.format(attributes[var][reg][bas])
                        else:
                            temp = '{:7.5f}'.format(attributes[var][reg][bas])
                        strout += '{:>10}'.format(temp)
                    f.write(strout+'\n')
                f.write('#-----------------------------------------------\n')
            f.close()
            if verbose: print("... written basis coefficient file to disk "+outfile)
        return


    def coefficients_to_cards(self, base_units = True):
        """evaluates bases coefficients at prescribed depth levels to get a card deck file

        base_units: convert from native units to base units in constants
        """
        import pint_pandas
        pint_pandas.PintType.ureg = constants.ureg

        # make an array of radii based on the first paramaterization that is read
        param_indx = 0
        radii = np.array([],dtype=float)
        for region in self.metadata['parameterization'][param_indx]:
            if region not in ['num_regions','filename','description']:
                top_temp = self.metadata['parameterization'][param_indx][region]['top_radius']
                bot_temp = self.metadata['parameterization'][param_indx][region]['bottom_radius']
                lvl_temp = self.metadata['parameterization'][param_indx][region]['levels']
                rad_temp = np.linspace(top_temp,bot_temp,num=lvl_temp)
                radii=np.append(radii,rad_temp)
        radii.sort() # sort from center to surface
        #names=['rho','vpv','vsv','qkappa','qmu','vph','vsh','eta']
        #use units from the elas/anelas file
        #names = tools.convert2nparray(self.metadata['parameters'])
        #units = tools.convert2nparray(self.metadata['units'])
        names=['rho','vpv','vsv','qkappa','qmu','vph','vsh','eta']
        units =['kg/m^3','m/s','m/s','dimensionless','dimensionless','m/s','m/s','dimensionless']


        fields=list(zip(names,units))

        # loop over names and call evaluate_at_depth
        # Create data array for converted to Panda array with units
        PA_ = pint_pandas.PintArray; temp_dict = {}; temp_dict['radius'] = PA_(radii, dtype="pint[km]")
        #self.evaluate_at_depth(24.4,'vsv')
        for paraindx,param in enumerate(names):
            val_temp = self.evaluate_at_depth(self.metadata['norm_radius']-radii,param)
            # overwrite the repeated radii at bottom
            bottom_indx = np.where(np.ediff1d(radii)==0)[0]
            val_temp2 = self.evaluate_at_depth(self.metadata['norm_radius']-radii[bottom_indx],param,boundary='-')
            for indx,val in enumerate(val_temp2): val_temp[bottom_indx[indx]] = val
            temp_dict[param] = PA_(val_temp, dtype="pint["+units[paraindx]+"]")
            # convert from fractions to absolute parameter (qkappa, qmu)
            if '/' in param: # convert from fractions such as 1000/qmu to qmu
                frac = param.split('/')
                numerator = float(frac[0])
                param = frac[-1]
                val_temp = np.divide(numerator,val_temp,out=np.zeros_like(val_temp), where=val_temp!=0)
            # loop over names and check if there's /name; modify units if needed
                temp_dict[param] = PA_(val_temp, dtype="pint[1/"+units[paraindx]+"]")
        modelarr = pd.DataFrame(temp_dict)
        if base_units: # convert to base units
            for col in modelarr.columns: modelarr[col] = modelarr[col].pint.to_base_units()
        modelarr['depth'] = PA_((constants.R.magnitude - modelarr['radius'].pint.to(constants.R.units).values.quantity.magnitude).tolist(), dtype = constants.R.units)

        self._nlayers = len(modelarr['radius'])
        self.data = modelarr
        self._radius_max = max(self.data['radius'])

    def write_mineos_cards(self,file):
        import pint_pandas

        if self.data is None or self._nlayers == 0: raise ValueError('reference1D data arrays are not allocated')
        names=['radius','rho','vpv','vsv','qkappa','qmu','vph','vsh','eta']
        units =['m','kg/m^3','m/s','m/s','dimensionless','dimensionless','m/s','m/s','dimensionless']

        # check if the units are the same or conversion is needed and where
        convert_columns = []
        for indx,name in enumerate(names):
            if pint_pandas.PintType.ureg.parse_expression(units[indx]).units != self.data[name].pint.units: convert_columns.append(name)

        if 'discontinuities' not in self.metadata.keys(): self.get_discontinuity()
        disc = self.metadata['discontinuities']
        # first write the header
        printstr  =  [unicode(self._name+"\n")]
        printstr.append(unicode("1 %.1f 1 1\n" % (self.metadata['ref_period'])))
        printstr.append(unicode("  %d  %d  %d  %d  %d\n" % (self._nlayers,disc['itopic']+1,disc['itopoc']+1,disc['itopmantle']+1,disc['itopcrust']+1)))

        shape = self.data[names].shape
        output = np.zeros(shape)
        for ii in range(shape[0]):
            for jj in range(shape[1]):
                if units[jj] not in convert_columns:
                    output[ii,jj] = self.data[names].iloc[ii,jj].to(units[jj]).magnitude
                else:
                    output[ii,jj] = self.data[names].iloc[ii,jj].magnitude
        # write the string in the fortran format
        header_line  =  ff.FortranRecordWriter('f8.0,3f9.2,2f9.1,2f9.2,f9.5')
        for ii in range(shape[0]): printstr.append(unicode(header_line.write(output[ii,:])+'\n'))
        printstr[-1] = printstr[-1].split('\n')[0]
        # write the file
        f= open(file,"w")
        f.writelines(printstr)
        f.close()
        return

    def read_mineos_cards(self,file,header = 3):
        # Operations between PintArrays of different unit registry will not work.
        # We can change the unit registry that will be used in creating new
        # PintArrays to prevent this issue.
        import pint_pandas
        pint_pandas.PintType.ureg = constants.ureg

        names=['radius','rho','vpv','vsv','qkappa','qmu','vph','vsh','eta']
        units =['m','kg/m^3','m/s','m/s','dimensionless','dimensionless','m/s','m/s','dimensionless']
        fields=list(zip(names,units))
        #formats=[float for ii in range(len(fields))]
        # modelarr = np.genfromtxt(file,dtype=None,comments='#',skip_header=3,names=fields)
        modelarr = pd.read_csv(file,skiprows=header,comment='#',sep='\s+',names=fields)
        # read the punit units from last header
        modelarr_ = modelarr.pint.quantify(level=-1)

        # Create data array
        PA_ = pint_pandas.PintArray
        modelarr_['depth'] = constants.R - modelarr_['radius'].pint.to(constants.R.units)
        self.data = modelarr_
        self._radius_max = max(self.data['radius'])
        self.metadata['parameters'] = names[1:]
        self.metadata['units'] = units[1:]
        # Get the other metadata from the first 3 line header
        with open(file,'r') as f:
            head = [next(f).strip('\n') for x in range(header)]
        self._name = head[0].strip()
        self.metadata['ref_period'] = float(head[1].split()[1])
        self.metadata['norm_radius'] = constants.R.to('km').magnitude

        # No original metadata found
        self.metadata['attributes'] = None

        # store rest of the metadata
        self.metadata['description'] = 'Read from '+file
        self.metadata['filename'] = file
        self._nlayers = len(modelarr['radius'])

    def get_Love_elastic(self):
        '''
        Get the Love parameters and Voigt averaged elastic properties with depth

        A,C,N,L,F: anisotropy elastic Love parameters

        kappa: bulk modulus

        mu: shear modulus

        vphi: bulk sound velocity

        xi: shear anisotropy ratio

        phi: P anisotropy ratio

        lambda: A Lame constant = kappa-2/3mu

        Zs, Zp: S and P impedances
        '''
        if self.data is None or self._nlayers == 0: raise ValueError('reference1D data arrays are not allocated')

        # Add data fields
        self.data['a'] = self.data['rho']*self.data['vph']**2
        self.data['c'] = self.data['rho']*self.data['vpv']**2
        self.data['n'] = self.data['rho']*self.data['vsh']**2
        self.data['l'] = self.data['rho']*self.data['vsv']**2
        self.data['f'] = self.data['eta']*(self.data['a']-2.*self.data['l'])

        # equivalent isotropic
        self.data['kappa'] = (4.0*(self.data['a']+self.data['f']-self.data['n'])+self.data['c'])/9.
        self.data['mu'] = (self.data['a']+self.data['c']-2.*self.data['f']+5.*self.data['n']+6.*self.data['l'])/15.
        self.data['vp'] = ((self.data['kappa']+4.*self.data['mu']/3.)/self.data['rho']).pow(0.5)
        self.data['vs'] = (self.data['mu']/self.data['rho']).pow(0.5)
        self.data['vphi'] = (self.data['kappa']/self.data['rho']).pow(0.5)
        self.data['lambda'] = self.data['kappa']-2./3.*self.data['mu']

        # anisotropy
        self.data['xi'] = (self.data['vsh'].div(self.data['vsv'])).pow(2)
        self.data['phi'] = (self.data['vpv'].div(self.data['vph'])).pow(2)

        self.data['as'] = ((self.data['vsh']-self.data['vsv']).div((self.data['vsh']+self.data['vsv'])/2.).replace(np.nan, 0)).pint.to('percent')
        self.data['ap'] = ((self.data['vph']-self.data['vpv']).div((self.data['vsh']+self.data['vsv'])/2.).replace(np.nan, 0)).pint.to('percent')

        # impedance contrasts
        self.data['Zp'] = self.data['vp']*self.data['rho']
        self.data['Zs'] = self.data['vs']*self.data['rho']

        # Add metadata
        for field in ['a','c','n','l','f','vp','vs','vphi','xi','phi','Zp', 'Zs','kappa','mu','as','ap']:
            self.metadata['parameters'].append(field)
            self.metadata['units'].append(str(self.data[field].pint.units))

    def get_mineralogical(self):
        '''
        Get the Love parameters and Voigt averaged elastic properties with depth

        gravity: gavity at each depth

        Brunt-Vaisala Frequency: Used for Bullen's parameter

        Bullen: Bullen's parameter

        pressure: pressure at each depth

        poisson: Poisson's ratio
        '''
        if self.data is None or self._nlayers == 0:
            warnings.warn('reference1D data arrays are not allocated. Trying to apply coefficients_to_cards from within get_mineralogical')
            self.coefficients_to_cards()

        if not set(['vs', 'vp']).issubset(self.data.keys()) :
            warnings.warn('Derived elastic parameters not allocated. Trying to apply get_Love_elastic from within get_mineralogical')
            self.get_Love_elastic()

        # Operations between PintArrays of different unit registry will not work.
        # We can change the unit registry that will be used in creating new
        # PintArrays to prevent this issue.
        import pint_pandas
        pint_pandas.PintType.ureg = constants.ureg

        if constants.planetpreferred == 'Earth':
            file = tools.get_filedir()+'/'+self._name+'.'+str(uuid.uuid4())
            # write a temporary cards file
            self.write_mineos_cards(file)

            layers = self._nlayers
            grav,vaisala,bullen,pressure = getbullen(file,layers,constants.omega.to_base_units().magnitude,constants.G.to_base_units().magnitude)

            # Add data fields
            PA_ = pint_pandas.PintArray
            self.data['gravity'] = PA_(grav, dtype="pint[m/s^2]")
            self.data['Brunt-Vaisala'] = PA_(vaisala, dtype="pint[Hz]")
            self.data['Bullen'] = PA_(bullen, dtype="pint[dimensionless]")
            self.data['pressure'] = PA_(pressure, dtype="pint[Pa]")

            # Add metadata
            for field in ['gravity','Brunt-Vaisala','Bullen','pressure']:
                self.metadata['parameters'].append(field)
                self.metadata['units'].append(str(self.data[field].pint.units))

            # Poisson ratio
            vsbyvp1d = self.data['vs'].div(self.data['vp'])
            self.data['poisson2'] = np.divide((1.-(2.*vsbyvp1d*vsbyvp1d)),(2.-(2.*vsbyvp1d*vsbyvp1d)))
            self.metadata['parameters'].append('poisson')
            self.metadata['units'].append('dimensionless')

            # delete a file
            with contextlib.suppress(FileNotFoundError): os.remove(file)

        else:
            warnings.warn(' mineralogical parameters not evaluated for '+constants.planetpreferred)

    def if_discontinuity(self,depth_in_km):
        """
        Returns whether a depth is a discontinuity.
        """
        depth_in_km = tools.convert2nparray(depth_in_km)
        disc_array = self.metadata['discontinuities']['delta']['depth'].pint.to('km').values.quantity.magnitude
        output = np.zeros_like(depth_in_km,dtype=bool)
        for idep,depth in enumerate(depth_in_km):
            output[idep] = np.any(np.isclose(disc_array,depth))
        return output if len(output)>1 else output[0]

    def get_discontinuity(self):
        '''
        Get values, average values and contrasts at discontinuities

        Output:
        ------

        Returns a structure self.metadata['disc'] that has three arrays:

        delta: containing absolute difference in parameters between smaller/larger radii

        average: containing absolute average parameters between smaller/larger radii

        contrasts: containing contrast in parameters (in %)
        '''
        if self.data is None or self._nlayers == 0: raise ValueError('reference1D data arrays are not allocated')

        disc_depths = [item.magnitude for item, count in Counter(self.data['depth']).items() if count > 1]
        disc = {}
# Create a named array for discontinuities

        for field in ['delta','average','contrast']: disc[field] = 0. * self.data.copy().drop(range(len(np.unique(disc_depths)),len(self.data)))
        # convert units to percent in contrast
        for param in self.data.columns:
            if param != 'radius' and param != 'depth':
                disc['contrast'][param] = (0.*disc['contrast'][param]/disc['contrast'][param][0]).pint.to('percent')

        for icount,depth in enumerate(disc_depths):
            depth_ = constants.ureg.Quantity(depth,self.data['depth'].pint.units)
            sel = self.data[self.data['depth']==depth_]
            for field in sel:
                if field == 'radius' or field == 'depth':
                    disc['delta'][field][icount] = sel[field].iat[0]
                    disc['average'][field][icount] = sel[field].iat[0]
                    disc['contrast'][field][icount] = sel[field].iat[0]
                else:
                    disc['delta'][field][icount] = sel[field].iat[0]-sel[field].iat[1]
                    disc['average'][field][icount] = 0.5*(sel[field].iat[0]+sel[field].iat[1])
                    ## contrasts need to be in %
                    disc['contrast'][field][icount] = (abs(disc['delta'][field][icount]) / disc['average'][field][icount]).to('percent')

        #---- try to find discontinuities
        discradii = disc['delta']['radius'].pint.to('km').values.quantity.magnitude
        if 'vp' not in self.data.keys() or 'vs' not in self.data.keys(): self.get_Love_elastic()
        vp = self.data['vp'].pint.to('km/s').values.quantity.magnitude
        vs = self.data['vs'].pint.to('km/s').values.quantity.magnitude
        radii = self.data['radius'].pint.to('km').values.quantity.magnitude

        discfind = disc['delta']['radius'][np.abs(1221.5-discradii)<25.].pint.to('km').values.quantity.magnitude
        if len(discfind) <= 0: # not found
            warnings.warn(" itopic not found")
        elif len(discfind) > 1: raise ValueError('get_discontinuity: multiple values within discontinuity limits')
        else:
            disc['itopic'] = np.where(radii==discfind[0])[0][0]

        discfind = disc['delta']['radius'][np.abs(3480.0-discradii)<25.].pint.to('km').values.quantity.magnitude
        if len(discfind) <= 0: # not found
            warnings.warn(" itopoc not found")
        elif len(discfind) > 1:
            raise ValueError('get_discontinuity: multiple values within discontinuity limits')
        else:
            disc['itopoc'] = np.where(radii == discfind[0])[0][0]

        ###   Top of crust
        discfind = np.where(np.logical_and(vp < 7.5,vs > 0.))[0]
        if len(discfind) > 0: disc['itopcrust'] = max(discfind)
        #discfind = disc['delta']['radius'][np.abs(6368.0-disc['delta']['radius']/1000.)<0.1]
#         if len(discfind) <= 0: # not found
#             warnings.warn(" itopcrust not found")
#         elif len(discfind) > 1:
#             raise ValueError('get_discontinuity: multiple values within discontinuity limits')
#         else:
            #disc['itopcrust'] = np.where(self.data['radius']==discfind[0])[0][1]

        itopmantle = min(np.where(vp < 7.5)[0])
        if itopmantle >0: disc['itopmantle'] = itopmantle-1
        self.metadata['discontinuities'] = disc

    def get_custom_parameter(self,parameters):
        '''
        Get the arrays of custom parameters defined in various Earth models
        '''
        import pint_pandas
        if self.data is not None and self._nlayers > 0:
            PA_ = pint_pandas.PintArray
            # convert to array for ease of looping
            if isinstance(parameters,string_types): parameters = np.array([parameters])
            if isinstance(parameters,list): parameters = np.array(parameters)
            for ii in np.arange(parameters.size):
                if parameters[ii] not in self.metadata['parameters']:
                    if 'as' in parameters[ii]:
                        # Add data fields
                        self.data[parameters[ii]] = self.data['as']
                        self.metadata['parameters'].append(parameters[ii])
                        self.metadata['units'].append('percent')
                    elif 'ap' in parameters[ii]:
                        self.data[parameters[ii]] = self.data['ap']
                        self.metadata['parameters'].append(parameters[ii])
                        self.metadata['units'].append('percent')
                    elif 'vs' in parameters[ii]:
                        self.data[parameters[ii]] = self.data['vs']
                        self.metadata['parameters'].append(parameters[ii])
                        self.metadata['units'].append('m/s')
                    elif 'vp' in parameters[ii]:
                        self.data[parameters[ii]] = self.data['vp']
                        self.metadata['parameters'].append(parameters[ii])
                        self.metadata['units'].append('m/s')
                    elif 'rho' in parameters[ii]:
                        self.data[parameters[ii]] = self.data['rho']
                        self.metadata['parameters'].append(parameters[ii])
                        self.metadata['units'].append('kg/m^3')
                    else:
                        raise NotImplementedError('parameter ',parameters[ii],' is not currently implemented in reference1D.get_custom_parameter')
        else:
            raise ValueError('reference1D object is not allocated')

    def evaluate_at_depth(self,depth_in_km,parameter='vsh',boundary='+',interpolation='linear'):
        '''
        Get the values of a parameter at a given depth

        boundary: + for value at larger radius at a discontinuity
        '''
        # need to update so it can deal with vectors
        if boundary not in ['+','-']: raise ValueError('boundary needs to be - or +')
        depth_in_km = tools.convert2nparray(depth_in_km)
        values = np.zeros_like(depth_in_km,dtype=float)
        parameters = tools.convert2nparray(self.metadata['parameters'])
        matching = [(indx,s) for indx,s in enumerate(parameters) if parameter in s]
        if len(matching) == 0 : raise KeyError('parameter '+parameter+' cannot be evaluated from reference1d instance')
        unit = self.metadata['units'][matching[0][0]]

        # check for 1000/qmu or 1000/qkappa
        divide1000=False
        if (matching[0][1] != parameter) and matching[0][1].startswith('1000/') and (parameter in ['qkappa','qmu']): divide1000 = True

        uniqueregions = {}
        target_region = np.empty_like(depth_in_km,dtype="U40"); target_region[:] = ''

        # detailed information about the native parameterization which went into the
        # inversion is available
        if self.metadata['attributes'] is not None:
        # check if norm_radius is within reasonable range
            if not 0.98*constants.R.to('km').magnitude <= self.metadata['norm_radius'] <= 1.02*constants.R.to('km').magnitude :
                raise ValueError('Normalizing radius not compatible')
            radius_in_km = self.metadata['norm_radius'] - depth_in_km
            param_indx = self.metadata['attributes'][matching[0][1]]['param_index']
            # finding target region in depth
            for region in self.metadata['parameterization'][param_indx]:
                if region not in ['num_regions','filename','description']:
                    # difference with top and bottom radii
                    difftop = self.metadata['parameterization'][param_indx][region]['top_radius'] - radius_in_km
                    diffbot = self.metadata['parameterization'][param_indx][region]['bottom_radius']-radius_in_km
                    dep_flag = difftop*diffbot

                    # within a region if dep_flag is neative
                    flag_array = (dep_flag < 0)
                    # if it is 0 then choose the flag based on boundary
                    findzeroindx = np.where(dep_flag==0.)[0]
                    if len(findzeroindx) != 0:
                        for indx in findzeroindx: #depth corresponds to the bottom of a region
                            if depth_in_km[indx] == 0: # surface of the solid earth
                                flag_array[indx] = True
                            else:
                                if diffbot[indx] == 0 and boundary == '+': flag_array[indx] = True
                                if difftop[indx] == 0 and boundary == '-': flag_array[indx] = True

                    for idx,flag in enumerate(flag_array):
                        if flag:
                            target_region[idx] = region
                            # store unique regions
                            if region not in [*uniqueregions]:
                                uniqueregions[region] = {'radius_range':
                                [self.metadata['parameterization'][param_indx][region]['bottom_radius'],
                                self.metadata['parameterization'][param_indx][region]['top_radius']],
                                'types': [*self.metadata['attributes'][matching[0][1]][region]]}
            if np.any(target_region == ''): raise ValueError('target regions not found')
            # create arguments for bases evaluations
            rnorm = self.metadata['norm_radius']
            #loop over all regions that are relevant to these depths
            for region  in [*uniqueregions]:
                # find all depth queries that are in this region
                indx = np.where(target_region==region)[0]
                radial_splines = False # assume that no splines are included
                choices = ['TOP', 'BOTTOM', 'CONSTANT', 'LINEAR', 'QUADRATIC', 'CUBIC']
                if not np.all([(key in choices or key.startswith('DEGREE')) for key in uniqueregions[region]['types']]): radial_splines = True
                if not radial_splines:
                    # evaluate value of the polynomial coefficients and derivative at each depth
                    vercof,dvercof = tools.bases.eval_polynomial(radius_in_km[indx],
                            uniqueregions[region]['radius_range'] ,rnorm,
                            uniqueregions[region]['types'])
                else:
                    spline_config = self.metadata['parameterization'][param_indx][region]['polynomial'].split(':')
                    if spline_config[0].strip() != 'SPLINEPNTS':
                        raise AssertionError('Polynomial type should be SPLINEPNTS in ' + region)
                    nsplines = int(spline_config[-1])
                    vercof1,dvercof1 = tools.bases.eval_splrem(radius_in_km[indx], uniqueregions[region]['radius_range'], nsplines)
                    # in case there's only one depth, make sure the shape of vercof1 fit vercof
                    vercof1 = vercof1.reshape([len(indx),nsplines])
                    dvercof1 = dvercof1.reshape([len(indx),nsplines])
                    # select the relevant splines
                    splindx = []; nonspltype = []; isspl = np.zeros(len(uniqueregions[region]['types']),dtype='bool')
                    for typekey,spltype in enumerate(uniqueregions[region]['types']):
                        if spltype.startswith('SPLINE'):
                            splindx.append(int(spltype.split()[-1])-1)
                            isspl[typekey] = True
                        else:
                            nonspltype.append(spltype)
                    if np.all(isspl):
                        vercof = vercof1; dvercof = dvercof1
                    else:
                        # evaluate other non-spline bases
                        # doesnt seem correct
                        vercof2,dvercof2 = tools.bases.eval_polynomial(radius_in_km[indx],
                                uniqueregions[region]['radius_range'] ,rnorm,nonspltype)
                        # combine polynomials and splines in the original order
                        vercof = np.zeros([len(indx),len(uniqueregions[region]['types'])]);
                        dvercof = np.zeros([len(indx),len(uniqueregions[region]['types'])]);
                        vercof[:,isspl] = vercof1[:,splindx]; dvercof[:,isspl] = dvercof1[:,splindx]
                        vercof[:,~isspl] = vercof2; dvercof[:,~isspl] = dvercof2
                # build up the coefficient array
                coef = []
                for key in uniqueregions[region]['types']:
                    coef.append(self.metadata['attributes'][matching[0][1]][region][key])
                temp = np.dot(vercof,np.array([coef]).T)
                for key, val in enumerate(indx):
                    if temp.ndim == 1:
                        values[val] = temp[key]
                    else:
                        values[val] = temp[key][0]
                    if divide1000 and values[val]!= 0.: values[val]=1000./values[val]

        # If detailed metadata regarding the basis parameterization is not available
        # interpolated based on the card deck file
        else:
            if self.data is not None:
                if parameter in self.metadata['parameters']:
                    modelval = self.data[parameter].values.quantity.magnitude
                    depth_array = constants.R.to('km').magnitude - self.data['radius'].pint.to('km').values.quantity.magnitude # in km
                    values = np.zeros_like(depth_in_km)
                    # Only select the region to interpolate
                    disc_depth = self.metadata['discontinuities']['delta']['depth'].pint.to('km').values.quantity.magnitude
                    disc_depth.sort()
                    # append a depth to get the index of discontinuity below
                    disc_depth = np.append(disc_depth,constants.R.to('km').magnitude)
                    disc_depth = np.append(0.,disc_depth)

                    for idep,depth in enumerate(depth_in_km):
                        # is it a discontinity
                        if self.if_discontinuity(depth):
                            findindices = np.nonzero(np.isclose(depth_array,depth))[0]
                            if boundary == '+': # larger radius, smaller depth
                                values[idep] = modelval[findindices[1]]
                            elif boundary == '-':
                                values[idep] = modelval[findindices[0]]
                        else: # not a boundary
                            # this is the bottom of the region under consideration
                            indxdisc = np.where(disc_depth-depth>=0.)[0][0]
                            if indxdisc == 0: indxdisc += 1 # surface is queried, so use the region below for interpolation
                            end=np.where(np.isclose(depth_array-disc_depth[indxdisc-1],0))[0][0]
                            if indxdisc != 1:
                                # leave the first and last since these are other regions
                                start=np.where(np.isclose(depth_array-disc_depth[indxdisc],0))[0][1]
                                indxdeps = np.arange(start,end+1)
                            else:
                                start=0
                                indxdeps = np.arange(start,end+1)
                            values[idep] = griddata(depth_array[indxdeps],modelval[indxdeps],depth,method=interpolation)

                else:
                    raise ValueError('parameter '+parameter+' not defined in array')
            else:
                raise ValueError('data in reference1D object is not allocated')
        return values*constants.ureg(unit)

    def apply_adams_williamson(self, OCO=True, ICO=False, rhotopOC=9.90344, jumprhoIC=0.6, step_km=1.):

        if ICO and not OCO: raise ValueERROR('Cannot apply Adams Wiliamson relationship for density to inner core without first applying to inner core')
        if self.data is None or self._nlayers == 0:
            warnings.warn('reference1D data arrays are not allocated. Trying to apply coefficients_to_cards from within apply_adams_williamson')
            self.coefficients_to_cards()
        if not set(['gravity', 'kappa']).issubset(self.data.keys()) :
            warnings.warn('Derived elastic parameters not allocated. Trying to apply get_mineralogical from within apply_adams_williamson')
            self.get_mineralogical()

        itopic = self.metadata['discontinuities']['itopic']
        itopoc = self.metadata['discontinuities']['itopoc']
        rr = self.data['radius'].pint.to('km').values.quantity.magnitude
        dr=self.data['radius'].pint.to('km').diff(periods=-1).values.quantity.magnitude
        g = self.data['gravity'].pint.to('m/s^2').values.quantity.magnitude
        rho = self.data['rho'].pint.to('g/cc').values.quantity.magnitude
        phi = (self.data['kappa']/self.data['rho']).pint.to('meter ** 2 / second ** 2').values.quantity.magnitude
        rnorm = self.metadata['norm_radius']

        logicOC = (rr<=rr[itopoc]) & \
                (rr >= rr[itopic]) & \
                (self.data['vsh'].values.quantity.magnitude == 0)
        logicIC = (rr <= rr[itopic]) & \
                (self.data['vsh'].values.quantity.magnitude != 0)

        if OCO:
            parameter='rho'
            region='outer core'

            attributes = self.metadata['attributes'][parameter][region]
            types = list(attributes.keys())
            param_indx = self.metadata['attributes'][parameter]['param_index']
            meta_region = self.metadata['parameterization'][param_indx][region]
            radius_range = [meta_region['bottom_radius'],meta_region['top_radius']]

            # # prescribe density at top of outer core
            rho_list=[rhotopOC]
            rhotopOC = constants.ureg(str(rhotopOC)+' g / cc')
            self.data['rho'].iloc[itopoc] = rhotopOC

            radii_list = np.linspace(rr[itopoc],rr[itopic],int((rr[itopoc]-rr[itopic])/step_km+1))
            vercof,dvercof = tools.bases.eval_polynomial(radii_list, radius_range, rnorm, types)

            tck_rho = splrep(rr[logicOC], rho[logicOC])
            tck_g = splrep(rr[logicOC], g[logicOC])
            tck_phi = splrep(rr[logicOC], phi[logicOC])
            for radius_query in radii_list[:-1]:
                rho_query=splev(radius_query, tck_rho)
                g_query=splev(radius_query, tck_g)
                phi_query=splev(radius_query, tck_phi)
                drhodr=-rho_query*g_query/phi_query
                rho_list.append(rho_list[-1]-drhodr*step_km*1000)

            # find best fitting coeffients
            ATA = np.dot(vercof.T, vercof)
            solution = np.dot(np.linalg.inv(ATA),np.dot(vercof.T, rho_list))

            for ii,key in enumerate(types):
                try:
                    attributes[key] = np.round(solution[ii][0],5)
                except:
                    attributes[key] = np.round(solution[ii],5)

        if ICO:
            parameter='rho'
            region='inner core'

            attributes = self.metadata['attributes'][parameter][region]
            types = list(attributes.keys())
            param_indx = self.metadata['attributes'][parameter]['param_index']
            meta_region = self.metadata['parameterization'][param_indx][region]
            radius_range = [meta_region['bottom_radius'],meta_region['top_radius']]

            rhotopIC = self.data['rho'].iloc[itopic+1]+constants.ureg(str(jumprhoIC)+' g / cc')
            # first calculate the values at every step
            rho_list=[rhotopIC.to('g/cc').magnitude]
            self.data['rho'].iloc[itopoc] = rhotopIC

            radii_list = np.linspace(rr[itopic],0,int(rr[itopic]/step_km+1))
            vercof,dvercof = tools.bases.eval_polynomial(radii_list, radius_range, rnorm, types)

            tck_rho = splrep(rr[logicIC], rho[logicIC])
            tck_g = splrep(rr[logicIC], g[logicIC])
            tck_phi = splrep(rr[logicIC], phi[logicIC])
            for radius_query in radii_list[:-1]:
                rho_query=splev(radius_query, tck_rho)
                g_query=splev(radius_query, tck_g)
                phi_query=splev(radius_query, tck_phi)
                drhodr=-rho_query*g_query/phi_query
                rho_list.append(rho_list[-1]-drhodr*step_km*1000)

            # find best fitting coeffients
            ATA = np.dot(vercof.T, vercof)
            solution = np.dot(np.linalg.inv(ATA),np.dot(vercof.T, rho_list))

            for ii,key in enumerate(types):
                try:
                    attributes[key] = np.round(solution[ii][0],5)
                except:
                    attributes[key] = np.round(solution[ii],5)

        # redo come calculations
        self.coefficients_to_cards()
        self.derive()

    def to_mineoscards(self,directory='.',fmt='cards'):
        '''
        Writes a model file that is compatible with MINEOS.
        '''
        parameters = ['radius','rho','vpv','vsv','qkappa','qmu','vph','vsh','eta']
        units =['m','kg/m^3','m/s','m/s','dimensionless','dimensionless','m/s','m/s','dimensionless']
        if self.data is not None and self._nlayers > 0:
            model_name = self._name
            ntotlev = self._nlayers
            itopic = self.metadata['discontinuities']['itopic']+1
            itopoc = self.metadata['discontinuities']['itopoc']+1
            itopmantle = self.metadata['discontinuities']['itopmantle']+1
            itopcrust = self.metadata['discontinuities']['itopcrust']+1

            outfile = directory+'/'+model_name+'.'+fmt
            f = open(outfile,'w')
            f.write(model_name+'\n')
            f.write('1 1. 1 1\n')
            line = ff.FortranRecordWriter('(5I5)')
            f.write(line.write([ntotlev,itopic,itopoc,itopmantle,itopcrust])+u'\n')
            line = ff.FortranRecordWriter('(f8.0,3f9.2,2f9.1,2f9.2,f9.5)')

            write = self.data[parameters]
            for i in range(self._nlayers):
                eachline = []
                for indx,param in enumerate(parameters):
                    eachline.append(write.iloc[i][param].to(units[indx]).magnitude)
                f.write(line.write(eachline)+u'\n')
            f.close()
            print('... written mineos file '+outfile)
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
            model_name = self._name
            outfile = directory+'/'+model_name+'.'+fmt
            f = open(outfile,'w')
            f.write('{} - P\n'.format(model_name))
            f.write('{} - S\n'.format(model_name))

            for i in range(self._nlayers):
                d_idx = self._nlayers-i-1
                f.write('{:2.4f}   {:2.4f}   {:2.4f}    {:2.4f}\n'.format(
                   (self._radius_max - self.data['radius'][d_idx]).to('km').magnitude,
                   self.data['vp'][d_idx].to('km/s').magnitude,
                   self.data['vs'][d_idx].to('km/s').magnitude,
                   self.data['rho'][d_idx].to('g/cm^3').magnitude))
            f.close()
            print('... written TauP file '+outfile)
        else:
            raise ValueError('reference1D object is not allocated')

        if self.data['vp'].iloc[-1] == 0 or self.data['vs'].iloc[-1] == 0:
            raise Warning('zero velocity layer detected at surface ...\n \
                      TauP raytracing may not work')

    def to_axisem(self,directory='.',anelastic=True,anisotropic=True,fmt='bm'):
        '''
         Write 1D model to be used as an external model in axisem
        '''
        if self.data is not None and self._nlayers > 0:
            model_name = self._name
            outfile = directory+'/'+model_name+'.'+fmt
            f = open(outfile,'w')
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

            keys = ['radius','rho','vpv','vsv','qkappa','qmu','vph','vsh','eta']
            values =['m','kg/m^3','m/s','m/s','dimensionless','dimensionless','m/s','m/s','dimensionless']
            units = dict(zip(keys, values))

            for i in range(self._nlayers):

                f.write('{}    {}    {}    {}    {}    {}    {}    {}    {}\n'.format(
                self.data['radius'][::-1][i].to(units['radius']).magnitude,
                self.data['rho'][::-1][i].to(units['rho']).magnitude,
                self.data['vpv'][::-1][i].to(units['vpv']).magnitude,
                self.data['vsv'][::-1][i].to(units['vsv']).magnitude,
                self.data['qkappa'][::-1][i].to(units['qkappa']).magnitude,
                self.data['qmu'][::-1][i].to(units['qmu']).magnitude,
                self.data['vph'][::-1][i].to(units['vph']).magnitude,
                self.data['vsh'][::-1][i].to(units['vsh']).magnitude,
                self.data['eta'][::-1][i].to(units['eta']).magnitude) )

                if i < self._nlayers-1 and self.data['radius'][::-1][i] == self.data['radius'][::-1][i+1]:
                    depth_here = (self._radius_max - self.data['radius'][::-1][i]) / 1000.0
                    n_discon += 1
                    f.write('#    Discontinuity {}, depth {:6.2f} km\n'.format(n_discon,depth_here))

            f.close()
            print('... written AXISEM file '+outfile)

        else:
            raise ValueError('reference1D object is not allocated')

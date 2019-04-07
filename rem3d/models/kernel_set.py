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
from scipy import sparse

####################### IMPORT REM3D LIBRARIES  #######################################
from .lateral_basis import lateral_basis
from .radial_basis import radial_basis
from .common import radial_attributes
from .. import tools
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
        self.initialize(dict)     
        self.extract_lateral(dict)
        self.extract_radial(dict)
        
    def initialize(self,dict,required = ['nmodkern','ivarkern','desckern','ncoefhor','ncoefcum','nhorpar','ihorpar','ityphpar','typehpar','numvar','varstr'],optional = ['forward_modeling','scaling']):
        for var in required:
            try:
                self.metadata[var] = dict[var]
            except:
                raise KeyError('required field '+var+' not found for kernel_set')
        for var in optional:
            try:
                self.metadata[var] = dict[var]
            except:
                self.metadata[var] = None
        
    def extract_lateral(self,dict):
        lateral=[]
        for ihor in np.arange(self.metadata['nhorpar']):
            type = self.metadata['typehpar'][ihor]
            metadata = {}
            metadata['ncoefhor']=dict['ncoefhor'][ihor]
            if 'SPHERICAL HARMONICS' in type:
                metadata['lmaxhor'] = dict['lmaxhor'][ihor]
            elif 'PIXELS' in type:
                for field in ['xsipix','xlapix','xlopix']:
                    metadata[field] = np.array(dict[field][ihor], order = 'F')
            elif 'SPHERICAL SPLINES' in type:
                for field in ['ixlspl','xlaspl','xlospl','xraspl']:
                    metadata[field] = np.array(dict[field][ihor], order = 'F')
            else:
                raise NotImplementedError(type+' has not been implemented in kernel_set.extract_lateral')
            lateral.append(lateral_basis(name='HPAR'+str(ihor+1), type = type, metadata=metadata))
        self.data['lateral_basis']=np.array(lateral)
        
    def extract_radial(self,dict):
        radial={}
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        for variable in dict['varstr']: #loop over all variables, grabbing
            radial[variable]=[]
            
            ivarfind =np.where(self.metadata['varstr']==variable)[0]
            assert(len(ivarfind) == 1),'only one parameter can be selected in eval_kernel_set'
            findrad = np.array([(ii, dict['desckern'][ii]) for ii in np.arange(len(dict['ivarkern'])) if ivarfind[0]+1 == self.metadata['ivarkern'][ii]],dtype=dt)

            metadata = {};found = False
            types = np.unique([findrad['kernel'][ii].split(',')[-2].strip() for ii in np.arange(len(findrad))])
            assert(len(types) == 1),'only one type is allowed'
            
            for jj in np.arange(len(findrad)):
                radker = findrad['kernel'][jj]
                metadata = {}; 
                found = False
                if 'variable splines' in radker or 'vbspl' in radker:
                    found = True
                    metadata['knots'] = [float(findrad['kernel'][ii].split(',')[-1].split('km')[0]) for ii in np.arange(len(findrad))]
                    radial[variable].append(radial_basis(name=radker, type = 'variable splines', metadata=metadata))
                    
                elif 'delta' in radker or 'dirac delta' in radker:
                    found = True
                    metadata['info'] = radker.split(',')[-1]
                    radial[variable].append(radial_basis(name=radker, type = 'dirac delta', metadata=metadata))
                elif 'boxcar' in radker or 'constant' in radker:
                    found = True                    
                    metadata['depthtop'] = [float(findrad['kernel'][ii].split(',')[-1].split('-')[0]) for ii in np.arange(len(findrad))]
                    metadata['depthbottom'] = [float(findrad['kernel'][ii].split(',')[-1].split('-')[1].split('km')[0]) for ii in np.arange(len(findrad))]
                    radial[variable].append(radial_basis(name=radker, type = 'boxcar', metadata=metadata))
                if not found: raise ValueError('information not found for '+radker)
        self.data['radial_basis']=radial
        
    def getprojection(self,latitude,longitude,depth_in_km,parameter='(SH+SV)*0.5'):
        
        # all the parameter options
        stringsum = self.metadata['varstr'][0]
        for stri in self.metadata['varstr'][1:]: stringsum= stringsum+', '+stri

        # select the radial kernels for this parameter
        dt = np.dtype([('index', np.int), ('kernel', np.unicode_,50)])
        ivarfind =np.where(self.metadata['varstr']==parameter)[0]
        assert(len(ivarfind) == 1),'only one parameter can be selected in eval_kernel_set among: '+stringsum+'. Only '+str(len(ivarfind))+' found for parameter '+parameter
        findrad = np.array([(ii, self.metadata['desckern'][ii]) for ii in np.arange(len(self.metadata['ivarkern'])) if ivarfind[0]+1 == self.metadata['ivarkern'][ii]],dtype=dt)

        # select corresponding lateral bases
        lateral_basis = self.data['lateral_basis']
        try:
            lateral_select = lateral_basis[self.metadata['ihorpar']-1][findrad['index']]
        except:
            raise ValueError('ihorpar needs to be defined for a kernel set. The HPAR for each radial kernel')
 
        #make sure only one variable is selected based on parameter input
        variables = np.unique(self.metadata['varstr'][self.metadata['ivarkern']-1][findrad['index']])
        assert(len(variables) == 1),'only one parameter, not '+str(len(variables))+', can be selected in eval_kernel_set from: '+stringsum        # select radial bases for this variable
        radial_select = self.data['radial_basis'][variables[0]]
        
        #initialize a projection matrix
        proj = sparse.csr_matrix((1,self.metadata['ncoefcum'][-1]))
        # loop over all radial kernels that belong to this parameter and add up
        for ii in np.arange(len(radial_select)):
            # start and end of indices to write to
            indend = self.metadata['ncoefcum'][findrad['index'][ii]]-1
            if findrad['index'][ii] == 0:
                indstart = 0
            else:
                indstart = self.metadata['ncoefcum'][findrad['index'][ii]-1]
            vercof, dvercof = radial_select[ii].eval_radial(depth_in_km)
            # convert to numpy arrays
            vercof = tools.convert2nparray(vercof)
            if vercof[ii] != 0.:
                horcof = lateral_select[ii].eval_lateral(latitude,longitude)
                proj=proj+sparse.csr_matrix( (horcof.data*vercof[ii],horcof.indices+indstart,horcof.indptr), shape=(1,self.metadata['ncoefcum'][-1]))
        return proj
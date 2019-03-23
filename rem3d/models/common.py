#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format."""

#####################  IMPORT STANDARD MODULES   ######################################   
# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

import sys,os
import timeit
import glob # pattern matching for finding files
import numpy as np #for numerical analysis
import pdb    #for the debugger pdb.set_trace()
import fortranformat as ff #reading/writing fortran formatted text
from configobj import ConfigObj
import xarray as xr
from io import StringIO
from copy import deepcopy
import ntpath
import warnings
import pandas as pd
import pint # For SI units
ureg = pint.UnitRegistry()

####################### IMPORT REM3D LIBRARIES  #######################################
from .. import tools 
from .reference1D import reference1D
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
    
def writeepixfile(filename,epixarr,headers=['#BASIS:PIX','#FORMAT:50']):
    """Write .epix file format from a named array.

    Parameters
    ----------

    filename : Name of the file containing four columns
              (latitude, longitude, pixel_size, value)

    """
    #combine headers
    header=''
    for hh in headers: header=header+'\n'+hh
    currentdir=os.getcwd()
    try:
        np.savetxt(filename, epixarr, fmt='%8.3f %8.3f %8.3f  %+12.7e',header=header,comments='')
    except :
        raise ValueError("File (",filename,") cannot be written in the current directory - ",currentdir)

    return 

def radial_attributes(desckern):
    """
    Takes a list of kernel descriptions and returns the attributes to be used in 
    radial_basis.py
    """    
    
def read3dmodelfile(modelfile):
    """
    Reads a standard 3D model file. maxkern is the maximum number of radial kernels
    and maxcoeff is the maximum number of corresponding lateral basis functions.
    resolution and realization are the indices for the resolution level
    and the realization from a model ensemble (usually 0 if a single file)
    """    
    if (not os.path.isfile(modelfile)): raise IOError("Filename ("+modelfile+") does not exist")

    #defaults
    desckern=[];ihorpar=[]
    refmodel = None; kerstr = None; crust = None; null_model= None
    interpolant = None; cite = None; shortcite = None; scaling = None
    with open(modelfile) as f: lines = f.readlines()
    ii=0;   model3d = {}
    while ii < len(lines):
        line=lines[ii]; ii=ii+1
        if line.startswith("REFERENCE MODEL:"): refmodel=line[16:].rstrip('\n').strip(' ')
        if line.startswith("NAME:"): model_name=line[5:].rstrip('\n').strip(' ')
        if line.startswith("KERNEL SET:"): kerstr=line[12:].rstrip('\n').strip(' ')
        if line.startswith("NULL MODEL:"): 
            null_model=line[12:].rstrip('\n').strip(' ')
            if 'None' in null_model or 'none' in null_model: null_model=None
        if line.startswith("CITE:"): cite=line[6:].rstrip('\n').strip(' ')
        if line.startswith("SHORTCITE:"): shortcite=line[11:].rstrip('\n').strip(' ')
        if line.startswith("INTERPOLANT:"): interpolant=line[13:].rstrip('\n').strip(' ')
        if line.startswith("CRUST:"): crust=line[7:].rstrip('\n').strip(' ')
        if line.startswith("SCALING:"): scaling=line[9:].rstrip('\n').strip(' ')
        if line.startswith("RADIAL STRUCTURE KERNELS:"): nmodkern = int(line[26:].rstrip('\n'))
        
        # search for parmaterization
        foundsplines = False; foundharmonics = False; foundpixels = False              
        if line.startswith("DESC"): 
            idummy=int(line[4:line.index(':')])
            substr=line[line.index(':')+1:len(line.rstrip('\n'))]
            ifst,ilst=tools.firstnonspaceindex(substr)
            desckern.append(substr[ifst:ilst])
        if line.startswith("HORIZONTAL PARAMETERIZATIONS:"): 
            nhorpar = int(line[29:].rstrip('\n'))
            typehpar=np.zeros(nhorpar, dtype='U40')
            ityphpar=np.zeros(nhorpar, dtype=np.int)
            ncoefhor=np.zeros(nhorpar, dtype=np.int)
            hsplfile=np.zeros(nhorpar, dtype='U40')
            coef=[[] for i in range(nmodkern)]
        if line.startswith("HPAR"):     
            idummy=int(line[4:line.index(':')])
            substr=line[line.index(':')+1:len(line.rstrip('\n'))]
            ifst,ilst=tools.firstnonspaceindex(substr)
            if substr[ifst:ifst+20] == 'SPHERICAL HARMONICS,':
                # initialize
                lmaxhor=np.zeros(nhorpar, dtype=np.int)
                
                ityphpar[idummy-1]=1
                typehpar[idummy-1]='SPHERICAL HARMONICS'
                lmax = int(substr[21:].rstrip('\n'))
                lmaxhor[idummy-1]=lmax
                ncoefhor[idummy-1]=(lmax+1)**2
                foundharmonics = True
            elif substr[ifst:ifst+18] == 'SPHERICAL SPLINES,': 
                ifst1=ifst+18
                ifst=len(substr)
                ilst=len(substr)
                while substr[ifst-1:ifst] != ',': ifst=ifst-1
                ncoef=int(substr[ifst+1:ilst].rstrip('\n'))
                substr=substr[ifst1:ifst-1]
                ifst1,ilst=tools.firstnonspaceindex(substr)
                hsplfile[idummy-1]=substr[ifst1:ilst]
                ityphpar[idummy-1]=2
                typehpar[idummy-1]='SPHERICAL SPLINES'
                ncoefhor[idummy-1]=ncoef
                
                # initialize if found first time
                if not foundsplines:
                    ixlspl=[[] for i in range(nhorpar)]
                    xlaspl=[[] for i in range(nhorpar)]
                    xlospl=[[] for i in range(nhorpar)]
                    xraspl=[[] for i in range(nhorpar)]
                
                # specific variables
                for jj in range(ncoef):
                    arr=lines[ii].rstrip('\n').split(); ii=ii+1
                    ixlspl[idummy-1].append(int(arr[0]))
                    xlaspl[idummy-1].append(float(arr[1]))
                    xlospl[idummy-1].append(float(arr[2]))
                    xraspl[idummy-1].append(float(arr[3]))
                foundsplines = True        
            elif substr[ifst:ifst+7] == 'PIXELS,':  
                ifst1=ifst+7
                ifst=len(substr)
                ilst=len(substr)
                while substr[ifst-1:ifst] != ',': ifst=ifst-1
                ncoef=int(substr[ifst+1:ilst].rstrip('\n'))
                substr=substr[ifst1:ifst-1]
                ifst1,ilst=tools.firstnonspaceindex(substr)
                hsplfile[idummy-1]=substr[ifst1:ilst]
                ityphpar[idummy-1]=3
                typehpar[idummy-1]='PIXELS'
                ncoefhor[idummy-1]=ncoef
                
                # initialize
                if not foundsplines:
                    xlapix=[[] for i in range(nhorpar)]
                    xlopix=[[] for i in range(nhorpar)]
                    xsipix=[[] for i in range(nhorpar)]
                
                # specific variables
                for jj in range(ncoef):
                    arr=lines[ii].rstrip('\n').split(); ii=ii+1
                    xlopix[idummy-1].append(float(arr[0]))
                    xlapix[idummy-1].append(float(arr[1]))
                    xsipix[idummy-1].append(float(arr[2]))
                foundpixels = True              
            else:
                raise ValueError('Undefined parameterization type - '+substr[ifst:ilst])
        if line.startswith("STRU"):
            idummy=int(line[4:line.index(':')])
            ihor=int(line[line.index(':')+1:].rstrip('\n'))
            ihorpar.append(ihor)
            ncoef=ncoefhor[ihor-1]
            for jj in range(int(ncoef/6)):
                arr=lines[ii].rstrip('\n').split(); ii=ii+1
                #coef[jj*6:(jj+1)*6,idummy-1]=[float(i) for i in arr]
                for i in arr: coef[idummy-1].append(float(i))
            remain = ncoef % 6    
            if remain > 0: 
                arr=lines[ii].rstrip('\n').split(); ii=ii+1
                for val in arr: coef[idummy-1].append(float(val))
    # Store the variables
    numvar=0; varstr=np.zeros(nmodkern, dtype='U40')
    ivarkern=np.zeros(nmodkern)
    for ii in np.arange(nmodkern):
        string=desckern[ii]
        #pdb.set_trace()
        if numvar == 0:
            varstr[0] = string[:string.index(',')]
            ivarkern[ii]=1
            numvar=numvar+1
        else:
            for kk in np.arange(numvar):
                if varstr[kk] == string[:string.index(',')]:
                    ivarkern[ii]=kk+1
        
        if ivarkern[ii] == 0:
            numvar=numvar+1
            varstr[numvar-1] = string[:string.index(',')]
            ivarkern[ii]=numvar

    # Save the relevant portions
    desckern = np.array(desckern)
    ihorpar = np.array(ihorpar)
    varstr = varstr[:numvar]
    ncoefcum = np.cumsum([ncoefhor[ihor-1] for ihor in ihorpar])
    
    # Store in a dictionary
    metadata = {}
    metadata['refmodel']=refmodel; metadata['kerstr']=kerstr;metadata['nmodkern']=nmodkern
    metadata['desckern']=desckern; metadata['nhorpar']=nhorpar;metadata['hsplfile']=hsplfile
    metadata['ityphpar']=ityphpar; metadata['typehpar']=typehpar; 
    metadata['ncoefhor']=ncoefhor; metadata['ihorpar']=ihorpar
    metadata['ivarkern']=ivarkern; metadata['numvar']=numvar
    metadata['varstr']=varstr; metadata['ncoefcum']=ncoefcum
    metadata['crust']=crust; metadata['null_model']=null_model
    metadata['interpolant']=interpolant; metadata['cite']=cite
    metadata['shortcite']=shortcite; metadata['scaling']=scaling

    if 'SPHERICAL HARMONICS' in typehpar:
        metadata['lmaxhor']=lmaxhor
    if 'PIXELS' in typehpar:
        metadata['xsipix']=np.array(xsipix); metadata['xlapix']=np.array(xlapix)
        metadata['xlopix']=np.array(xlopix)  
    if 'SPHERICAL SPLINES' in typehpar:
        metadata['ixlspl']=np.array(ixlspl); metadata['xlaspl']=np.array(xlaspl)
        metadata['xlospl']=np.array(xlospl); metadata['xraspl']=np.array(xraspl)
        
    model3d['data']={}
    model3d['data']['coef']=pd.DataFrame(coef)
    try:
        model3d['data']['name']=model_name
    except:
        model3d['data']['name']=ntpath.basename(modelfile)
    model3d['metadata']=metadata
    return model3d
    
def readensemble(filename):
    """Read .npz file containing a collection of models.

    Parameters
    ----------

    filename : An ensemble file

    """


    return 

def readprojections(type='radial'):
    """Read .npz file containing a collection of models.

    Parameters
    ----------

    filename : An file containing all projection matrices

    """


    return 
    
def epix2xarray(model_dir='.',setup_file='setup.cfg',output_dir='.',n_hpar=1,write_zeros=True,buffer=True):
    """
    Input parameters:
    ----------------
  
    model_dir: path to directory containing epix layer files
    
    output_dir: path to output directory
        
    setup_file: setup file containing metadata for the model
           
    n_hpar: number of unique horizontal parameterizations
    
    buffer: write to buffer instead of file for the intermediate step of the ascii file

    """

    start_time  =  timeit.default_timer()  
    
    cfg_file = model_dir+'/'+setup_file

    if not os.path.isfile(cfg_file):
        raise IOError('No configuration file found.'\
	                 'Model directory must contain '+setup_file)
    else:
        parser = ConfigObj(cfg_file)
    model_name = parser['metadata']['name']
    kernel_set = '{}'.format(parser['metadata']['kerstr'])
      
    if buffer:
        print('... writing ASCII buffer')
    else:
        print('... writing ASCII file')
    asciibuffer = epix2ascii(model_dir=model_dir,setup_file=setup_file,output_dir=output_dir,n_hpar=n_hpar,write_zeros=write_zeros,buffer=buffer)
    elapsed  =  timeit.default_timer() - start_time
    print("........ evaluations took "+str(elapsed)+" s")
    if buffer:
        print('... read to ASCII buffer. evaluations took '+str(elapsed)+' s')
    else:
        print('... written ASCII file '+asciibuffer+'. evaluations took '+str(elapsed)+' s')
    ncfile = output_dir+'/{}.{}.rem3d.nc4'.format(model_name,kernel_set)
    print('... writing netcdf file '+ncfile)
    ds = ascii2xarray(asciibuffer,outfile=ncfile,setup_file=setup_file)
    
    return ds
    
def epix2ascii(model_dir='.',setup_file='setup.cfg',output_dir='.',n_hpar=1,write_zeros=True, checks=True,buffer=False, onlyheaders=False):
    '''
    write a rem3d formatted ascii file from a directory containing epix files 

    Input parameters:
    ----------------
  
    epix_dir: path to directory containing epix layer files
    
    setup_file: setup file containing metadata for the model

    output_file: name of rem3d format output file
        
    n_hpar:number of horizontal parameterizations (currently only handles 
           models with 1 horizontal parameterization, and with a constant pixel width)
           
    checks: if True, checks if the metadata in setup_file is consistent with epix files

    buffer: write to buffer instead of file for the intermediate step of the ascii file
    
    onlyheaders: only write headers, not the coefficients
    '''
    cfg_file = model_dir+'/'+setup_file
    ref_dict = {} #dictionary containing reference values

    if not os.path.isfile(cfg_file):
        raise IOError('No configuration file found.'\
	                 'Model directory must contain '+setup_file)
    else:
        parser = ConfigObj(cfg_file)

    model_name = parser['metadata']['name']
    ref_model = parser['metadata']['refmodel']
    epix_folder = parser['metadata']['folder']
    interpolant = parser['metadata']['interpolant']
    cite = parser['metadata']['cite']
    crust = parser['metadata']['crust']
    forward_modeling = parser['metadata']['forward_modeling']
    scaling = parser['metadata']['scaling']
    kernel_set = '{}'.format(parser['metadata']['kerstr'])

    #write header
    outfile = output_dir+'/{}.{}.rem3d.ascii'.format(model_name,kernel_set)
    if buffer:
        # Writing to a buffer
        f_out = StringIO()
    else:
        f_out = open(outfile,'w')
    f_out.write(u'NAME: {}\n'.format(model_name))
    f_out.write(u'REFERENCE MODEL: {} \n'.format(ref_model))
    f_out.write(u'KERNEL SET: {}\n'.format(kernel_set))
    f_out.write(u'INTERPOLANT: {}\n'.format(interpolant))
    f_out.write(u'CITE: {}\n'.format(cite))
    if crust is not 'None': f_out.write(u'CRUST: {}\n'.format(crust))
    f_out.write(u'FORWARD MODELING: {}\n'.format(forward_modeling))
    f_out.write(u'SCALING: {}\n'.format(scaling))

    #find the number radial kernels
    epix_lengths = []; string = []
    icount = 0 
    for parameter in parser['parameters']:
        mod_type = parser['parameters'][parameter]['type']
        par_folder = parser['parameters'][parameter]['folder']
        description = parser['parameters'][parameter]['description']
        icount += 1
        string.append(str(icount)+'. '+description+' ('+parameter+')')

        if mod_type == 'heterogeneity':
            epix_files = glob.glob(model_dir+'/'+epix_folder+'/'+par_folder+'/*.epix')
        elif mod_type == 'topography':
            epix_files = glob.glob(model_dir+'/'+epix_folder+'/'+par_folder+'/*'+parameter+'.epix')

        if len(epix_files)>0:
            epix_lengths.append(len(epix_files))
        else:
            raise IOError('no files found for parameter '+parameter+' of type '+mod_type)
    print('... read '+str(np.sum(epix_lengths))+' radial structure kernels of '+str(len(string))+' variables: \n'+'\n'.join(string))
    f_out.write(u'RADIAL STRUCTURE KERNELS: {}\n'.format(np.sum(epix_lengths)))

    n_hpar = 1 #default is a single horizontal parameterization for all parameters

    stru_indx = []
    stru_list = []
    lats = []
    lons = []
    pxs = []

    k = 1
    for i, parameter in enumerate(parser['parameters']):
        ref_dict[parameter] = {}
        ref_dict[parameter]['ifremav'] = []
        ref_dict[parameter]['refvalue'] = []
        ref_dict[parameter]['average'] = []

        mod_type = parser['parameters'][parameter]['type']
        par_folder = parser['parameters'][parameter]['folder']

        if mod_type == 'heterogeneity':
            epix_files = glob.glob(model_dir+'/'+epix_folder+'/'+par_folder+'/*.epix')
            epix_files.sort(key=tools.alphanum_key)
            ref_dict[parameter]['depth_in_km'] = []
        elif mod_type == 'topography':
            #topo_folder = parser['parameters'][parameter]['folder']
            epix_files = glob.glob(model_dir+'/'+epix_folder+'/'+par_folder+'/*'+parameter+'.epix')
            depth = parser['parameters'][parameter]['depth']*ureg(parser['parameters'][parameter]['unit'])
            depth.ito('km')
            ref_dict[parameter]['depth_in_km'] = float(depth.magnitude)
        else:
            raise ValueError('model type not recognized... should be either "heterogeneity" or "topography"')

        #write descriptions of radial kernels
        for j, epix_file in enumerate(epix_files):

            #read header and store basic metadata
            metadata = {}
            head = []
            with open(epix_file) as f:
                for line in f:
                    if line.startswith('#'):
                        head.append(line)
                    for field in ['DEPTH_IN_KM','AVERAGE','IFREMAV','REFVALUE','REFMODEL','UNIT','WHAT','FORMAT','BASIS']:
                        if field in line: metadata[field] = line.split(':')[1].split('\n')[0].lstrip().rstrip()
                                    
            # conduct checks
            if checks:
                assert (parser['parameters'][parameter]['unit'].lower()==metadata['UNIT'].lower())," in file "+epix_file                                

                                               
                if parameter.lower()!=metadata['WHAT'].lower()  or parser['parameters'][parameter]['description'].lower() == metadata['WHAT'].lower():
                    warnings.warn("parameter or its description !=metadata['WHAT'] in file "+epix_file)
                if parser['metadata']['refmodel'].lower()!=metadata['REFMODEL'].lower():
                    warnings.warn("parser['parameters']['refmodel']!=metadata['REFMODEL'] in file "+epix_file)
                assert (metadata['FORMAT']=='50')," in file "+epix_file
                assert (metadata['BASIS'].lower()=='PIX'.lower())," in file "+epix_file

            # defaults if field not available in the epix file
            try:
                ref_dict[parameter]['ifremav'].append(np.float(metadata['IFREMAV']))
            except:
                ref_dict[parameter]['ifremav'].append(0.)
            try:
                ref_dict[parameter]['refvalue'].append(np.float(metadata['REFVALUE']))
            except:
                ref_dict[parameter]['refvalue'].append(-999.0)
            try:
                ref_dict[parameter]['average'].append(np.float(metadata['AVERAGE']))
            except:
                ref_dict[parameter]['average'].append(0.)
            try:
                ref_dict[parameter]['refmodel'] = metadata['REFMODEL']
            except:
                ref_dict[parameter]['refmodel'] = 'None'
            
            if mod_type == 'heterogeneity':
                for line in head:
                    if 'DEPTH_RANGE' in line: depth_range = line.split(':')[1].split('\n')[0] 
                f_out.write(u'DESC  {:3.0f}: {}, boxcar, {} km\n'.format(k,parameter,depth_range))
                ref_dict[parameter]['depth_in_km'].append( np.float(metadata['DEPTH_IN_KM']))
            elif mod_type == 'topography':
                if checks: assert (float(parser['parameters'][parameter]['depth']) == float(metadata['REFVALUE']))," in file "+epix_file
                depth_ref = parser['parameters'][parameter]['depth']
                f_out.write(u'DESC  {:3.0f}: {}, delta, {} km\n'.format(k,parameter,depth_ref))
            
            # now read the data
            f = np.loadtxt(epix_file)

            if j == 0 and k == 1:
                par1 = f[:,0:2]
                lats.append(f[:,0])
                lons.append(f[:,1])
                pxs.append(f[:,2])
                stru_indx.append(n_hpar)
            else:
                par_new = f[:,0:2]
                lats_new = f[:,0]
                lons_new = f[:,1]
                pxs_new = f[:,2]

                #loop through all previous parameterizations
                for ii in range(len(lats)):

                    if not np.array_equal(par1,par_new):
                        n_hpar += 1
                        lats.append(f[:,0])
                        lons.append(f[:,1])
                        pxs.append(f[:,2])
                stru_indx.append(n_hpar)
            k += 1
            
            #enforce longitudes from 0 to 360 to be consistent with xarray
            assert(min(f[:,1]) >= 0.)," longitudes need to be [0,360] "+epix_file


    #write horizontal parameterization
    f_out.write(u'HORIZONTAL PARAMETERIZATIONS: {}\n'.format(len(lats)))
    for i in range(0,len(lats)):
        
        #check pixel widths
        if np.min(pxs[i]) != np.max(pxs[i]):
            raise ValueError('no inhomogeneous model parameterizations allowed yet')
        else:
            px_w = pxs[i][0]

        shape = (int(180.0/px_w),int(360.0/px_w))
        lats_arr = np.reshape(lats,shape,order='F')
        lons_arr = np.reshape(lons,shape,order='F')
        px_w_arr = np.reshape(pxs,shape,order='F')

        lats_ = lats_arr.flatten()
        lons_ = lons_arr.flatten()
        pxs_ = px_w_arr.flatten()

        f_out.write(u'HPAR   {}: PIXELS,  {:3.2f} X {:3.2f}, {}\n'.format(stru_indx[0],px_w,px_w,len(lats_)))

        assert (np.all(sorted(np.unique(lons_))==np.unique(lons_)))
        assert (np.all(sorted(np.unique(lats_))==np.unique(lats_)))
        for j in range(0,len(lats_)):
            lon_here = lons_[j]
            lat_here = lats_[j]
            px_here = pxs_[j]

            f_out.write(u'{:6.2f} {:6.2f} {:6.2f}\n'.format(lon_here,lat_here, px_here))
    
    
    if not onlyheaders:
        # write coefficients
        k = 1
        for i, parameter in enumerate(parser['parameters']):
            #epix_files = glob.glob(model_dir+'/'+epix_folder+'/'+parameter+'/*.epix')
            mod_type = parser['parameters'][parameter]['type']
            par_folder = parser['parameters'][parameter]['folder']

            if mod_type == 'heterogeneity':
                epix_files = glob.glob(model_dir+'/'+epix_folder+'/'+par_folder+'/*.epix')
                epix_files.sort(key=tools.alphanum_key)
            elif mod_type == 'topography':
                #topo_folder = parser['parameters'][parameter]['folder']
                epix_files = glob.glob(model_dir+'/'+epix_folder+'/'+par_folder+'/*'+parameter+'.epix')
            else:
                raise ValueError('model type not recognized... should be either "heterogeneity" or "topography"')

            # read the 1D model if any of the reference values are not defined
            ifread1D = np.any(np.array(ref_dict[parameter]['refvalue'])<0.)
            if ifread1D: 
                try: # try reading the 1D file in card format 
                    ref1d = reference1D(ref_dict[parameter]['refmodel'])
                    if mod_type == 'heterogeneity': ref1d.get_custom_parameter(parameter)
                except:
                    ifread1D = False
        
            #write model coefficients
            line = ff.FortranRecordWriter('(6E12.4)')
            for j, epix_file in enumerate(epix_files):
                f = np.loadtxt(epix_file)
                print('writing coefficients for layer ', k)
                coefs = f[:,3]
            
                #check if the reference value is negative. 
                # if so, make an instance of the 1D
                # model class to read from
                if ifread1D and ref_dict[parameter]['refvalue'][j] < 0: 
                    depth_in_km = ref_dict[parameter]['depth_in_km'][j]
                    ref_dict[parameter]['refvalue'][j] = ref1d.evaluate_at_depth(depth_in_km,parameter)
                    
                #check ifremav. if it's 1, add in average
                if ref_dict[parameter]['ifremav'][j] == 1:
                    coefs += refs_dict[parameter]['average'][j]
                    print('... adding average back to parameter '+parameter+' # '+str(j))

                coefs_arr = np.reshape(coefs,shape,order='F')
                coefs = coefs_arr.flatten()
                f_out.write(u'STRU  {:3.0f}:  {:1.0f}\n'.format(k,px_w))
                f_out.write(line.write(coefs)+u'\n')
                k += 1
    if buffer:
        f_out.seek(0)
        return f_out
    else:
        return outfile

def ascii2xarray(asciioutput,outfile=None,setup_file='setup.cfg',complevel=9, engine='netcdf4', writenc4 = True):
    '''
    write an xarrary dataset from a rem3d formatted ascii file

    Input parameters:
    ----------------
  
    ascii_file: path to rem3d format output file
        
    outfile: output netcdf file
    
    setup_file: setup file containing metadata for the model
    
    complevel, engine: options for compression in netcdf file
    
    writenc4: write a netcdf4 file, if True

    '''

    model_dict = {}
    
    # check for configuration file
    if not os.path.isfile(setup_file):
        raise IOError('No configuration file found.'\
                     'Model directory must contain '+setup_file)
    else:
        parser = ConfigObj(setup_file)

    try: #attempt buffer
        asciioutput.seek(0)
    except:
        if outfile == None: outfile = asciioutput+'.nc4'
        asciioutput = open(asciioutput,'r')
        
    #read header
    line = asciioutput.readline()
    while line:
        if 'NAME' in line:
            model_name = line.split('NAME:')[1].strip()
        if 'REFERENCE MODEL' in line:
            ref_model = line.split('REFERENCE MODEL:')[1].strip()
        if 'KERNEL SET' in line:
            krnl_set = line.split('KERNEL SET:')[1].strip()
        if 'RADIAL STRUCTURE KERNELS' in line:
            nrad_krnl = line.split('RADIAL STRUCTURE KERNELS:')[1].strip()
            nrad_krnl = int(nrad_krnl)
            break
        line = asciioutput.readline()
        
    # check that reference model is the same as parser
    assert(ref_model == parser['metadata']['refmodel']),ref_model+' the reference model in '+asciioutput+' is not the same as refmodel in '+setup_file
    assert(krnl_set == parser['metadata']['kerstr']) ,krnl_set+' the kernel string in '+asciioutput+' is not the same as kerstr in '+setup_file


    #read variables and parameterizations
    variables = []
    rpar_list = []
    hpar_list = []
    rpar = []
    variable_idxs = []

    npts_dep = nrad_krnl

    hpar_idx = 0
    rpar_idx = -1
    var_idx = 0

    i = 0
    line = asciioutput.readline()
    while line:
        if i < nrad_krnl:

            variable = line.strip().split()[2].split(',')[0]
            if variable not in variables:
                variables.append(variable)
                model_dict[variable] = {}
                model_dict[variable]['hpar_idx'] = None
                variable_idxs.append(var_idx)

                if len(rpar) > 0 and rpar not in rpar_list:
                    rpar_idx += 1
                    rpar_list.append(rpar)
                    model_dict[variables[var_idx]]['rpar_idx'] = rpar_idx 
                    var_idx += 1

                if len(rpar) > 0 and rpar in rpar_list:
                    rpar_list.append(rpar)
                    model_dict[variables[var_idx]]['rpar_idx'] = rpar_idx 
                    var_idx += 1
                    rpar = []
                
            try:
                rpar_start = float(line.strip().split(',')[-1].split('-')[0].strip('km'))
                rpar_end = float(line.strip().split(',')[-1].split('-')[1].strip('km'))
                rpar.append((rpar_start + rpar_end)/2.)
            except IndexError:
                model_dict[variable]['rpar_idx'] = None
            line = asciioutput.readline()
            i += 1

        if i == nrad_krnl: 
            #read number of horizontal parameterizations
            nhpar = int(line.strip().split()[-1])
            nhpar = line.strip().split()[-1]
            break
            
    # check that information on variables in ascii file exists in setup.cfg
    for var in variables: assert(var in parser['parameters'].keys()),var+' not found as shortname in '+setup_file

    for i in range(0,nhpar):

        lons = []
        lats = []
        pxwd = []

        line = asciioutput.readline()
        print(line.strip())
        hpar_type = line.strip().split()[2].split(',')[0]
        hpar_name = line.split(':')[1].strip()

        if hpar_name.lower().startswith('pixel'):
             pxw_lon = float(line.strip().split()[3].strip(','))
             pxw_lat = float(line.strip().split()[5].strip(','))
             nlines = int(360.0/pxw_lon) * int(180/pxw_lat)
             assert(nlines == float(line.strip().split()[6].strip(',')))
        else:
            raise ValueError('only PIXEL parameterizations enabled')

        for j in range(nlines):
            line = asciioutput.readline()
            lons.append(float(line.strip().split()[0]))
            lats.append(float(line.strip().split()[1]))
            pxwd.append(float(line.strip().split()[2]))

    hpar_list.append([lons,lats,pxwd])

    #read coefficients
    for variable in variables:
        stru_idx = model_dict[variable]['rpar_idx']
        model_dict[variable]['layers'] = {}

        if stru_idx is not None:
            n_stru = len(rpar_list[stru_idx])
        else:
            n_stru = 1

        for i in range(0,n_stru):
            layer_coefs = []

            #get structure info
            line = asciioutput.readline()
            nstru = int(line.strip().split()[1].split(':')[0])
            nhparam = int(line.strip().split()[2])
            npts = len(hpar_list[nhparam-1][0][:])
            nlines = int(npts/6)

            if model_dict[variable]['hpar_idx'] == None:
                model_dict[variable]['hpar_idx'] = nhparam-1

            for j in range(nlines):
                 line = asciioutput.readline()
                 for coef in line.strip().split():
                     layer_coefs.append(float(coef))
            # if not a multiple of 6, add remainder        
            remain = nlines % 6    
            if remain > 0: 
                line = asciioutput.readline()
                for coef in line.strip().split(): 
                    layer_coefs.append(float(coef))

            model_dict[variable]['layers'][i] = layer_coefs

    #open xarray dataset
    ds = xr.Dataset()

    #make DataArrays for each variable, and add to the dataset
    for variable in variables:
        hpar_idx = model_dict[variable]['hpar_idx'] 
        pxw = hpar_list[hpar_idx][2][0]
        print(variable,': PXW', pxw)

        #create dims arrays
        lon = np.arange((pxw/2.),360.,pxw)
        lat = np.arange(-90.+(pxw/2.),90,pxw)

        stru_idx = model_dict[variable]['rpar_idx']
        
        if stru_idx is not None:
            dep = rpar_list[stru_idx]
            data_array = xr.DataArray(np.zeros((len(dep),len(lat),len(lon))),
                                      dims = ['depth','latitude','longitude'],
                                      coords=[dep,lat,lon])
            for i,layer in enumerate(model_dict[variable]['layers']):
                data_array[i,:,:] = np.reshape(model_dict[variable]['layers'][layer],
                                    (len(lat),len(lon)),order='C')
        else:
            data_array = xr.DataArray(np.zeros((len(lat),len(lon))),
                                      dims = ['latitude','longitude'],
                                      coords = [lat,lon])
            data_array[:,:] = np.reshape(model_dict[variable]['layers'][0],
                                    (len(lat),len(lon)),order='C')
        #-------------------------------------------------------------------------
        #add reference values at each depth as metadata to the Data_Array
        #-------------------------------------------------------------------------
        av_attrs = {}
        for keys in parser['parameters'][variable].keys():
            av_attrs[keys] = parser['parameters'][variable][keys].decode('utf-8')
        # read the 1D model if any of the reference values are not defined
        av_attrs['refmodel'] = parser['metadata']['refmodel']
        
        ifread1D = True
        try: # try reading the 1D file in card format 
            ref1d = reference1D(av_attrs['refmodel'])
        except:
            ifread1D = False
        
        if len(data_array.shape) == 3: # if 3-D variable
            # get the variable values
            if ifread1D: ref1d.get_custom_parameter(variable)
            av_depth = deepcopy(data_array.depth.values)
            refvalue = []; avgvalue = []
            for depth in av_depth: 
                if ifread1D: 
                    refvalue.append(ref1d.evaluate_at_depth(depth,parameter=variable))
                else:
                    refvalue.append(-999.0)
                # select the appropriate map
                mapval = data_array.sel(depth=depth)
                # get the average
                globalav,area,percentarea = tools.MeanDataArray(mapval)
                avgvalue.append(globalav)    
            av_attrs['refvalue'] = np.array(refvalue)
            av_attrs['average'] = np.array(avgvalue)
            
        #add Data_Array object to Data_Set
        data_array.attrs = av_attrs
        ds[variable] = data_array

    #Add overall attributes
    attrs = {}
    for key in parser['metadata'].keys():
        attrs[key] = parser['metadata'][key].decode('utf-8')
    ds.attrs = attrs
 
    # write to netcdf
    comp = {'zlib': True, 'complevel': complevel}
    encoding = {var: comp for var in ds.data_vars}
    if writenc4 and outfile != None: ds.to_netcdf(outfile,engine=engine,encoding=encoding)
    
    return ds
        

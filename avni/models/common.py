#!/usr/bin/env python
"""This script/module contains routines that are used to Earth models"""

#####################  IMPORT STANDARD MODULES   ######################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int

import os
import timeit
import glob # pattern matching for finding files
import numpy as np #for numerical analysis
import fortranformat as ff #reading/writing fortran formatted text
from configobj import ConfigObj
import xarray as xr
from io import StringIO
from copy import deepcopy
import ntpath
import warnings
import pandas as pd
import struct
import traceback
import typing as tp

####################### IMPORT AVNI LIBRARIES  #######################################
from .. import tools
from .. import constants
from .reference1d import Reference1D
#######################################################################################

def readepixfile(filename: str) -> tp.Tuple[np.ndarray, dict, list]:
    """Read a local file in extended pixel (.epix) format.

    Parameters
    ----------
    filename : str
        Name of the file containing four columns:
        (`latitude`, `longitude`, `pixel_size`, `value`)

    Returns
    -------
    tp.Tuple[np.ndarray, dict, list]
        First element is an array containing (`latitude`, `longitude`, `pixel_size`, `value`).

        Second element are metadata from input fields if specified.

        Third element are all other comments except lines containing metadata.

    Raises
    ------
    IOError
        File not found in local directory

    Examples
    --------
    >>> epixarr,metadata,comments = readepixfile('file.epix')

    Notes
    -----
    This function assumes cell-centered values within a pixel of size `pixel_size`.
    In order to interpolate in arbitrary locations, this assumption needs to be used.

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # Read the file if available
    try:
        f = open(filename, 'r')
        epixarr=np.genfromtxt(filename, dtype=None,comments="#",names=['latitude','longitude','pixel_size','value'])
    except IOError:
        raise IOError("File (",filename,") cannot be read.")

    # strip header and store basic metadata
    metadata = {}
    comments = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                if ':' in line:
                    field = line.lstrip('#').split(':')[0].strip()
                    metadata[field] = line.lstrip('#').split(':')[1].strip()
                else:
                    comments.append(line.strip())

    return epixarr,metadata,comments

def writeepixfile(filename: str, epixarr: np.ndarray ,
                  metadata: dict = {'BASIS':'PIX','FORMAT':'50'},
                  comments: list = None) -> None:
    """Write named numpy array to extended pixel format (.epix) file.

    Parameters
    ----------
    filename : str
        Name of the file containing four columns:
        (`latitude`, `longitude`, `pixel_size`, `value`)
    epixarr : np.ndarray
        array containing (`latitude`, `longitude`, `pixel_size`, `value`)
    metadata : dict, optional
        metadata from input fields if specified, by default {'BASIS':'PIX','FORMAT':'50'}
    comments : list, optional
        all other comments except lines containing metadata, by default None

    Raises
    ------
    IOError
        File cannot be written in local directory

    Examples
    --------
    >>> writeepixfile('file.epix',epixarr,metadata,comments)

    Notes
    -----
    This function assumes cell-centered values within a pixel of size `pixel_size`.
    In order to interpolate in arbitrary locations, this assumption needs to be used.

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # combine headers
    header=''
    for key in sorted(metadata.keys()): header=header+'#'+key+':'+metadata[key]+'\n'
    if comments is not None:
        for comment in comments: header = header+comment+'\n'
    header = header[:-1] # to get rid of the last \n in header

    # try writing the file
    try:
        np.savetxt(filename, epixarr, fmt='%8.3f %8.3f %8.3f  %+12.7e',header=header,comments='')
    except :
        raise IOError("File (",filename,") cannot be written.")

    return

def read3dmodelfile(modelfile: str) -> dict:
    """Reads a standard 3D model file into a `model3d` dictionary containing `data`
    and `metadata`. The `data` field contains the basis coefficients of the
    parameterization employed in the 3D model file. The `metadata` fields contain
    various parameters relevant for evaluating the 3D model at various locations.

    Parameters
    ----------
    modelfile : str
        Text file describing the 3D model

    Returns
    -------
    dict
        A `model3d` dictionary containing `data` and `metadata`

    Raises
    ------
    IOError
        File not found in local directory
    ValueError
        Parameterization type has not been implemented yet.
        Current ones include spherical harmonics, spherical splines and pixels.

    Examples
    --------
    >>> model3d = read3dmodelfile('S362ANI+M')

    Notes
    -----
    This function has several fields in the `data` and `matadata` fields.
    The details depend on what kind of parameterization is employed.
    In general, `maxkern` is the maximum number of radial kernels
    and `maxcoeff` is the maximum number of corresponding lateral basis functions,
    `resolution` and `realization` are the indices for the resolution level
    and the realization from a model ensemble (usually 0 if it is a single file).

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # check if file exists
    if (not os.path.isfile(modelfile)): raise IOError("Filename ("+modelfile+") does not exist")

    #defaults
    desckern=[];ihorpar=[]
    refmodel = None; kerstr = None; crust = None; null_model= None
    interpolant = None; cite = None; shortcite = None

    # read the fields
    with open(modelfile) as f: lines = f.readlines()
    ii=0;   model3d = {}
    foundsplines = False; foundpixels = False; #foundharmonics = False
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
        if line.startswith("RADIAL STRUCTURE KERNELS:"): nmodkern = int(line[26:].rstrip('\n'))

        # search for parmaterization
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
                #foundharmonics = True
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
                if not foundpixels:
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
    numvar=0; varstr=np.zeros(nmodkern, dtype='U200')
    ivarkern=np.zeros(nmodkern,dtype=np.int)
    for ii in np.arange(nmodkern):
        string=desckern[ii]
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
    desckern = np.array(desckern, dtype='U200')
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
    metadata['shortcite']=shortcite

    # Get number of coefficients - lmax, number of splines or pixels
    if 'SPHERICAL HARMONICS' in typehpar:
        metadata['lmaxhor']=lmaxhor
    if 'PIXELS' in typehpar:
        metadata['xsipix']=np.array(xsipix); metadata['xlapix']=np.array(xlapix)
        metadata['xlopix']=np.array(xlopix)
    if 'SPHERICAL SPLINES' in typehpar:
        metadata['ixlspl']=np.array(ixlspl); metadata['xlaspl']=np.array(xlaspl)
        metadata['xlospl']=np.array(xlospl); metadata['xraspl']=np.array(xraspl)

    # Fill the basis coefficients as pandas DataFrame
    model3d['data']={}
    model3d['data']['coef']=pd.DataFrame(coef)
    try:
        model3d['data']['name']=model_name
    except:
        model3d['data']['name']=ntpath.basename(modelfile)
    model3d['metadata']=metadata

    return model3d

def epix2xarray(model_dir: str = '.', setup_file: str = 'setup.cfg',
                output_dir: str = '.', n_hpar: int = 1, buffer: bool = True) -> xr.Dataset:
    """Convert a set of extended pixel format (.epix) files in a directory
    to a netCDF4 file in AVNI format using xarray.

    Parameters
    ----------
    model_dir : str, optional
        path to directory containing epix layer files, by default '.'
    setup_file : str, optional
        setup file containing metadata for the model, by default 'setup.cfg'
    output_dir : str, optional
        path to output directory, by default '.'
    n_hpar : int, optional
        number of unique horizontal parameterizations, by default 1
    buffer : bool, optional
        write to buffer instead of file for the intermediate step of the ascii file, by default True

    Returns
    -------
    xr.Dataset
        An xarray Dataset that stores the values from a 3D model.

    Raises
    ------
    IOError
        File not found in local directory

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    start_time  =  timeit.default_timer()

    cfg_file = model_dir+'/'+setup_file

    if not os.path.isfile(cfg_file):
        raise IOError('No configuration file found.'\
	                 'Model directory must contain '+setup_file)
    else:
        parser = ConfigObj(cfg_file)

    model_name = parser['metadata']['name']
    try:
        kernel_set = '{}'.format(parser['metadata']['kerstr']).strip()
    except KeyError:
        kernel_set = 'NATIVE'

    if buffer:
        print('... writing ASCII buffer')
    else:
        print('... writing ASCII file')
    asciibuffer = epix2ascii(model_dir=model_dir,setup_file=setup_file,output_dir=output_dir,n_hpar=n_hpar,buffer=buffer)
    elapsed  =  timeit.default_timer() - start_time
    print("........ evaluations took "+str(elapsed)+" s")
    if buffer:
        print('... read to ASCII buffer. evaluations took '+str(elapsed)+' s')
    else:
        print('... written ASCII file '+asciibuffer+'. evaluations took '+str(elapsed)+' s')
    ncfile = output_dir+'/{}.{}.avni.nc4'.format(model_name,kernel_set)
    print('... writing netcdf file '+ncfile)
    ds = ascii2xarray(asciibuffer,model_dir=model_dir,outfile=ncfile,setup_file=setup_file)

    return ds

def checksetup(parser):
    """Check the arguments in the model setup file and populate the default options

    Parameters
    ----------
    parser : A ConfigObj object

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # Required arguments
    for key in ['name','cite']:
        if key not in parser['metadata'].keys():
            raise KeyError(key+' needs to be specified')
    for field in ['min','units','max','resolution']:
        for type in ['lat','lon']:
            key = 'geospatial_'+type+'_'+field
            if key in parser['metadata'].keys():
                parser['metadata'][key] = parser['metadata'][key] if field=='units' else float(parser['metadata'][key])
                if field == 'unit':
                    try:
                        constants.ureg(key)
                    except:
                        raise ValueError('need to define '+key+' in acceptable units (% not allowed, use percent')
            else:
                raise KeyError(key+' needs to be specified in setup file')
    for field in ['min','units','max']:
        key = 'geospatial_vertical_'+field
        if key in parser['metadata'].keys():
            parser['metadata'][key] = parser['metadata'][key] if field=='units' else float(parser['metadata'][key])
        else:
            raise KeyError(key+' needs to be specified')

    # Optional arguments
    try:
        value = parser['metadata']['geospatial_vertical_positive'].strip()
    except KeyError:
        parser['metadata']['geospatial_vertical_positive'] = 'down'
    if parser['metadata']['geospatial_vertical_positive'].lower() not in ['up','down']:
        raise KeyError('geospatial_vertical_positive type can only be up or down')
    try:
        value = parser['metadata']['interpolant'].strip()
        parser['metadata']['interpolant'] = None if value.lower() == 'none' else value
    except KeyError:
        parser['metadata']['interpolant'] = 'nearest'
    if parser['metadata']['interpolant'] != None:
        if parser['metadata']['interpolant'].lower() not in ['nearest','smooth']:
            raise KeyError('interpolant type can only be nearest or smooth for compatibility with KDtree queries')
    try:
        value = parser['metadata']['refmodel'].strip()
        parser['metadata']['refmodel'] = None if value.lower() == 'none' else value
    except KeyError:
        parser['metadata']['refmodel'] = None
    try:
        value = parser['metadata']['crust'].strip()
        parser['metadata']['crust'] = None if value.lower() == 'none' else value
    except KeyError:
        parser['metadata']['crust'] = None
    try:
        value = parser['metadata']['null_model'].strip()
        parser['metadata']['null_model'] = None if value.lower() == 'none' else value
    except KeyError:
        parser['metadata']['null_model'] = None
    try:
        value = '{}'.format(parser['metadata']['kerstr']).strip()
        parser['metadata']['kerstr'] = 'NATIVE' if value.lower() == 'none' else value
    except KeyError:
        parser['metadata']['kerstr'] = 'NATIVE'
    try:
        value = parser['metadata']['forward_modeling'].strip()
        parser['metadata']['forward_modeling'] = parser['parameters'].keys() if value.lower() == 'none' else value.split(',')
    except KeyError:
        parser['metadata']['forward_modeling'] = parser['parameters'].keys()

    # check units
    for var in parser['parameters'].keys():
        for key in ['unit','absolute_unit']:
            if key in parser['parameters'][var].keys():
                value =parser['parameters'][var][key]
                try:
                    constants.ureg(value)
                except:
                    raise ValueError('need to define '+key+' for variable '+var+' in acceptable units, not '+value+' (% NOT allowed, use percent).')
            else:
                if key == 'unit': raise ValueError('need to define '+key+' for variable '+var+' in acceptable units (% NOT allowed, use percent).')

def epix2ascii(model_dir: str = '.', setup_file: str = 'setup.cfg',
               output_dir: str = '.', n_hpar: int = 1,
               checks: bool = True, buffer: bool = True, onlyheaders: bool = False):
    """Write AVNI formatted ASCII file from a directory containing extended pixel format (.epix) files

    Parameters
    ----------
    model_dir : str, optional
        path to directory containing epix layer files, by default '.'
    setup_file : str, optional
        setup file containing metadata for the model, by default 'setup.cfg'
    output_dir : str, optional
        path to output directory, by default '.'
    n_hpar : int, optional
        number of horizontal parameterizations (currently only handles
        models with 1 horizontal parameterization, and with a constant pixel width), by default 1
    checks : bool, optional
        checks if the metadata in `setup_file` is consistent with epix files, by default True
    buffer : bool, optional
        write to buffer instead of file for the intermediate step of the ascii file, by default True
    onlyheaders : bool, optional
        write only headers, not the coefficients, by default False

    Returns
    -------
    buffer or file name
        The memory buffer containing the file output (buffer=True) or the output file name on disk.

    Raises
    ------
    IOError
        Some epix files have not been found in the directory
    AssertionError
        Checks for compatibility across epix files are not satisfied
    ValueError
        Invalid values are found for some metadata fields

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # Check the configuration file
    cfg_file = model_dir+'/'+setup_file
    ref_dict = {} #dictionary containing reference values
    if not os.path.isfile(cfg_file):
        raise IOError('No configuration file found.'\
	                 'Model directory must contain '+setup_file)
    else:
        parser = ConfigObj(cfg_file)

    # Required arguments
    checksetup(parser)
    model_name = parser['metadata']['name']
    cite = parser['metadata']['cite']
    epix_folder = parser['metadata']['folder']
    kernel_set = parser['metadata']['kerstr']

    # get the optional arguments
    interpolant = parser['metadata']['interpolant']
    ref_model = parser['metadata']['refmodel']
    crust = parser['metadata']['crust']
    null_model = parser['metadata']['null_model']
    forward_modeling = parser['metadata']['forward_modeling']

    #write header
    outfile = output_dir+'/{}.{}.avni.ascii'.format(model_name,kernel_set)
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
    if crust != 'None': f_out.write(u'CRUST: {}\n'.format(crust))
    if forward_modeling != 'None': f_out.write(u'FORWARD MODELING: {}\n'.format(forward_modeling))
    if null_model != 'None': f_out.write(u'NULL_MODEL: {}\n'.format(null_model))

    #find the number radial kernels
    epix_lengths = []; string = []
    icount = 0
    for parameter in parser['parameters']:
        mod_type = parser['parameters'][parameter]['type']
        par_folder = parser['parameters'][parameter]['folder']
        icount += 1
        try:
            description = parser['parameters'][parameter]['description']
            string.append(str(icount)+'. '+description+' ('+parameter+')')
        except KeyError:
            string.append(str(icount)+'. '+parameter)

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

    # Read through various parameters
    stru_indx = []
    #stru_list = []
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
            depth = parser['parameters'][parameter]['depth']*constants.ureg(parser['parameters'][parameter]['unit'])
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
                        if field in line: metadata[field] = line.split(':')[1].split('\n')[0].strip()

            # conduct checks
            if checks:
                if not  parser['parameters'][parameter]['unit'].lower()==metadata['UNIT'].lower(): raise AssertionError("UNIT incompatible in file "+epix_file)


                if parameter.lower()!=metadata['WHAT'].lower() :
                    warnings.warn("parameter !=metadata['WHAT'] in file "+epix_file)
                value = metadata['REFMODEL'].strip()
                refmodel_epix = None if value.lower() == 'none' else value
                if not ref_model == refmodel_epix:
                    warnings.warn("parser['parameters']['refmodel']!=metadata['REFMODEL'] in file "+epix_file)
                if not metadata['FORMAT']=='50': raise AssertionError("FORMAT incompatible in file "+epix_file)
                if not metadata['BASIS'].lower()=='PIX'.lower(): raise AssertionError("BASIS incompatible in file "+epix_file)

            # defaults if field not available in the epix file
            try:
                ref_dict[parameter]['ifremav'].append(float(metadata['IFREMAV']))
            except:
                ref_dict[parameter]['ifremav'].append(0.)
            try:
                ref_dict[parameter]['refvalue'].append(float(metadata['REFVALUE']))
            except:
                ref_dict[parameter]['refvalue'].append(-999.0)
            try:
                ref_dict[parameter]['average'].append(float(metadata['AVERAGE']))
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
                ref_dict[parameter]['depth_in_km'].append(float(metadata['DEPTH_IN_KM']))
            elif mod_type == 'topography':
                if checks:
                    if not float(parser['parameters'][parameter]['depth']) == float(metadata['REFVALUE']):
                        raise AssertionError("REFVALUE incompatible in file "+epix_file)
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
                        lats.append(lats_new)
                        lons.append(lons_new)
                        pxs.append(pxs_new)
                stru_indx.append(n_hpar)
            k += 1

            #enforce longitudes from 0 to 360 to be consistent with xarray
            if not min(f[:,1]) >= 0.: raise AssertionError("longitudes need to be [0,360] "+epix_file)

    #write horizontal parameterization
    f_out.write(u'HORIZONTAL PARAMETERIZATIONS: {}\n'.format(len(lats)))
    for i,_ in enumerate(lats):

        #check pixel widths
        if np.min(pxs[i]) != np.max(pxs[i]):
            raise ValueError('no inhomogeneous model parameterizations allowed yet')
        else:
            px_w = pxs[i][0]

        f_out.write(u'HPAR   {}: PIXELS,  {:5.3f} X {:5.3f}, {}\n'.format(stru_indx[0],px_w,px_w,len(lats[i])))

        if not np.all(sorted(np.unique(lons))==np.unique(lons)): raise AssertionError()
        if not np.all(sorted(np.unique(lats))==np.unique(lats)): raise AssertionError()
        for j in range(len(lats[i])):
            lon_here = lons[i][j]
            lat_here = lats[i][j]
            px_here = pxs[i][j]
            f_out.write(u'{:7.3f} {:7.3f} {:7.3f}\n'.format(lon_here,lat_here, px_here))

    if not onlyheaders:
        # read the 1D model if any of the reference values are not defined
        ifread1D = {} #stores whether the references values have been read from a file
        fileread = False
        for _, parameter in enumerate(parser['parameters']):
            mod_type = parser['parameters'][parameter]['type']
            ifread1D[parameter] = np.any(np.array(ref_dict[parameter]['refvalue'])<0.)
            if ifread1D[parameter]:
                # check if the refmodel file exists
                if not os.path.isfile(ref_dict[parameter]['refmodel']):
                    ifread1D[parameter] = False
                    print ('WARNING: Could not fill reference values for '+parameter+' as the 1D reference model file could not be found : '+ref_dict[parameter]['refmodel'])
            if ifread1D[parameter] :
                if not fileread:
                    try: # try reading the 1D file in card format
                        ref1d = Reference1D(ref_dict[parameter]['refmodel'])
                        fileread = True
                    except:
                        print ('WARNING: Could not fill reference values for '+parameter+' as the 1D reference model file could not be read as Reference1D instance : '+ref_dict[parameter]['refmodel'])
                        ifread1D[parameter] = False
            if mod_type == 'heterogeneity' and ifread1D[parameter]: ref1d.get_custom_parameter(parameter)

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

            #write model coefficients
            print('writing coefficients for '+parameter+' # '+str(k)+' - '+str(k+len(epix_files)-1)+' out of '+str(np.sum(epix_lengths))+' radial kernels/layers.')
            line = ff.FortranRecordWriter('(6E12.4)')
            for j in range(len(epix_files)):
                f = np.loadtxt(epix_files[j])
                coefs = f[:,3]

                #check if the reference value is negative.
                # if so, make an instance of the 1D
                # model class to read from
                if ifread1D[parameter] and ref_dict[parameter]['refvalue'][j] < 0:
                    depth_in_km = ref_dict[parameter]['depth_in_km'][j]
                    ref_dict[parameter]['refvalue'][j] = ref1d.evaluate_at_depth(depth_in_km,parameter)

                #check ifremav. if it's 1, add in average
                if ref_dict[parameter]['ifremav'][j] == 1:
                    coefs += ref_dict[parameter]['average'][j]
                    print('... adding average back to parameter '+parameter+' # '+str(j))
                f_out.write(u'STRU  {:3.0f}:  {:1.0f}\n'.format(k,px_w))
                f_out.write(line.write(coefs)+u'\n')
                k += 1

    # Write to buffer or file name
    if buffer:
        f_out.seek(0)
        return f_out
    else:
        return outfile

def ascii2xarray(asciioutput,outfile = None, model_dir: str = '.',
                 setup_file: str = 'setup.cfg', complevel: int = 9,
                 engine: str = 'netcdf4', writenc4: bool = False):
    """Create an xarray Dataset from AVNI formatted ASCII file

    Parameters
    ----------
    asciioutput : buffer or str
        A filename string or buffer that will be read
    outfile : str, optional
        Output file in netCDF4 format, by default None
    model_dir : str, optional
        _description_, by default '.'
    setup_file : str, optional
        setup file containing metadata for the model, by default 'setup.cfg'
    complevel : int, optional
        options for compression in netcdf file, by default 9
    engine : str, optional
        options in netcdf file, by default 'netcdf4'
    writenc4 : bool, optional
        write a netcdf4 file, by default False

    Returns
    -------
    ds: xarray Dataset
        Contains the model description in a new xarray Dataset

    Raises
    ------
    IOError
        Configuration file has not been found in the directory
    AssertionError
        Checks for compatibility across files are not satisfied
    ValueError
        Only pixels are allowed as parameterization for netCDF files
    warnings.warn
        Checks for number of pixels not satisfied

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    model_dict = {}
    cfg_file = model_dir+'/'+setup_file

    # check for configuration file
    if not os.path.isfile(cfg_file):
        raise IOError('No configuration file found.'\
                     'Model directory must contain '+setup_file)
    else:
        parser = ConfigObj(cfg_file)

    # Required arguments
    epix_folder = parser['metadata']['folder']
    checksetup(parser)

    try: #attempt buffer
        asciioutput.seek(0)
    except:
        if outfile == None and writenc4:
            try:
                outfile = asciioutput.split('.ascii')[0]+'.nc4'
            except:
                outfile = asciioutput+'.nc4'
        asciioutput = open(asciioutput,'r')

    #read header
    line = asciioutput.readline()
    while line:
        #if 'NAME' in line:
        #    model_name = line.split('NAME:')[1].strip()
        if 'REFERENCE MODEL' in line:
            value = line.split('REFERENCE MODEL:')[1].strip()
            ref_model = None if value.lower() == 'none' else value
        if 'KERNEL SET' in line:
            value = line.split('KERNEL SET:')[1].strip()
            krnl_set = None if value.lower() == 'none' else value
        if 'RADIAL STRUCTURE KERNELS' in line:
            nrad_krnl = line.split('RADIAL STRUCTURE KERNELS:')[1].strip()
            nrad_krnl = int(nrad_krnl)
            break
        line = asciioutput.readline()

    # check that reference model is the same as parser
    if not ref_model == parser['metadata']['refmodel']:
        raise AssertionError(ref_model+' the reference model in '+asciioutput+' is not the same as refmodel in '+setup_file)
    if not krnl_set == parser['metadata']['kerstr']: raise AssertionError(krnl_set+' the kernel string in '+asciioutput+' is not the same as kerstr in '+setup_file)

    #read variables and parameterizations
    variables = []
    rpar_list = []
    hpar_list = []
    rpar_starts = {}
    rpar_ends = {}
    rpar = {}
    variable_idxs = []
    hpar_idx = 0
    rpar_idx = 0
    var_idx = 0

    i = 0
    line = asciioutput.readline()
    while line:
        if i < nrad_krnl:
            variable = line.strip().split()[2].split(',')[0]
            if variable not in variables:
                variables.append(variable)
                rpar_starts[variable] = []
                rpar_ends[variable] = []
                rpar[variable] = []
                model_dict[variable] = {}
                model_dict[variable]['hpar_idx'] = None
                variable_idxs.append(var_idx)
                var_idx += 1

            # first try to find start and end of radial param
            try:
                rpar_start = float(line.strip().split(',')[-1].split('-')[0].strip('km'))
                rpar_end = float(line.strip().split(',')[-1].split('-')[1].strip('km'))
                rpar_starts[variable].append(rpar_start)
                rpar_ends[variable].append(rpar_end)
                rpar[variable].append((rpar_start + rpar_end)/2.)
            except IndexError:
                model_dict[variable]['rpar_idx'] = None
            line = asciioutput.readline()
            i += 1

        if i == nrad_krnl:
            #read number of horizontal parameterizations
            nhpar = int(line.strip().split()[-1])
            break

    # Now get rparindex
    for variable in variables:
        if len(rpar[variable]) != 0: # if it is an empty list like in discontinuity
            if len(rpar[variable]) > 1:
                if sorted(rpar[variable]) != rpar[variable]: raise AssertionError('depths not sorted',rpar[variable])
            if rpar[variable] not in rpar_list:
                rpar_list.append(rpar[variable])
                model_dict[variable]['rpar_idx'] = rpar_idx
                rpar_idx += 1
            else:
                model_dict[variable]['rpar_idx'] = rpar_list.index(rpar[variable])

    # check that information on variables in ascii file exists in setup.cfg
    for var in variables:
        if not var in parser['parameters'].keys(): raise AssertionError(var+' not found as shortname in '+setup_file)
        for indx in ['rpar_idx','hpar_idx']:
            if not indx in model_dict[var].keys(): raise AssertionError(var+' not read properly with index '+indx)

    for i in range(nhpar):
        lons = []
        lats = []
        pxwd = []

        line = asciioutput.readline()
        print(line.strip())
        #hpar_type = line.strip().split()[2].split(',')[0]
        hpar_name = line.split(':')[1].strip()

        if hpar_name.lower().startswith('pixel'):
             pxw_lon = float(line.strip().split()[3].strip(','))
             pxw_lat = float(line.strip().split()[5].strip(','))
             nlines_input = int(line.strip().split()[6].strip(','))
             nlines = int(360.0/pxw_lon) * int(180/pxw_lat)
             if not nlines == nlines_input: warnings.warn('number of pixels expected for '+str(pxw_lat)+'X'+str(pxw_lon)+' is '+str(nlines)+',  not '+str(nlines_input)+' as provided.')
        else:
            raise ValueError('only PIXEL parameterizations enabled')

        for j in range(nlines_input):
            line = asciioutput.readline()
            lons.append(float(line.strip().split()[0]))
            lats.append(float(line.strip().split()[1]))
            pxwd.append(float(line.strip().split()[2]))

        hpar_list.append([lons,lats,pxwd])

    #read coefficients
    num_coefficients = None
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
            #nstru = int(line.strip().split()[1].split(':')[0])
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

            # Check if the same number of coefficients are used
            if num_coefficients is None:
                num_coefficients = len(model_dict[variable]['layers'][i])
            else:
                if num_coefficients != len(model_dict[variable]['layers'][i]):
                    raise AssertionError('number of coefficients for variable '+variable+' and layer '+str(i)+' is not equal to that of another one '+str(num_coefficients)+'. All variables should use the same grid.')

    # check if we can  read 1D model
    ifread1D = False if parser['metadata']['refmodel'] is None else True
    if ifread1D:
        if os.path.isfile(parser['metadata']['refmodel']):
            try: # try reading the 1D file in card format
                ref1d = Reference1D(parser['metadata']['refmodel'])
                # get derived parameters
                ref1d.get_mineralogical()
                ifread1D = True
            except:
                ifread1D = False
                print ('WARNING: Could not fill some reference values as the 1D reference model file could not be read as Reference1D instance : '+parser['metadata']['refmodel'])
        else:
            ifread1D = False
            print ('WARNING: Could not fill some reference values as the 1D reference model file could not be found : '+parser['metadata']['refmodel'])

    # check that indices have been read properly
    for var in variables:
        if model_dict[var]['hpar_idx'] is None: raise AssertionError(var+' not read properly with hpar_idx ')

    # Get all depths
    alldepths = []; allstartdepths = []; allenddepths = []
    for rpar_temp in rpar_list: alldepths.extend(rpar_temp)
    alldepths = np.sort(np.unique(np.asarray(alldepths)))
    for variable in rpar_starts.keys():
        allstartdepths.extend(rpar_starts[variable])
        allenddepths.extend(rpar_ends[variable])
    allstartdepths = np.sort(np.unique(np.asarray(allstartdepths)))
    allenddepths = np.sort(np.unique(np.asarray(allenddepths)))

    #open xarray dataset
    ds = xr.Dataset()

    #make DataArrays for each variable, and add to the dataset
    area = None # calculate area the first time around
    for variable in variables:
        hpar_idx = model_dict[variable]['hpar_idx']

        # sort them by increaing lon if needed since reshape requires that
        sortbylon = np.all(hpar_list[hpar_idx][0] == sorted(hpar_list[hpar_idx][0]))

        # unique values for reshaping
        arr=pd.DataFrame(np.asarray(hpar_list[hpar_idx]).T,columns =['lon', 'lat', 'pxw'])
        if sortbylon:
            lon = np.unique(hpar_list[hpar_idx][0])
            lat = np.unique(hpar_list[hpar_idx][1])
            pxw = np.unique(hpar_list[hpar_idx][2])
        else:
            arr = arr.sort_values(by=['lon','lat'])
            lon = pd.unique(arr['lon'])
            lat = pd.unique(arr['lat'])
            pxw = pd.unique(arr['pxw'])

        if not len(pxw)==1: warnings.warn('more than 1 pixel size in variable '+variable)
        print(variable,': PXW', pxw)

        #create dims arrays
        stru_idx = model_dict[variable]['rpar_idx']

        # get the grid sizes stored
        pixel_array = xr.DataArray(np.zeros((len(lat),len(lon))),
                        dims = ['latitude','longitude'],
                        coords = [lat,lon])
        pixel_array[:,:] = np.reshape(arr['pxw'].values,
                            (len(lat),len(lon)),order='F')
        ds['pixel_width'] = pixel_array

        if stru_idx is not None:
            dep = rpar_list[stru_idx]
            # find depth indices
            dep_indx = np.searchsorted(alldepths,dep)
            data_array = xr.DataArray(np.zeros((len(alldepths),len(lat),len(lon))),
                                      dims = ['depth','latitude','longitude'],
                                      coords=[alldepths,lat,lon])
            for i,layer in enumerate(model_dict[variable]['layers']):
                # if sorting is needed before reshaping
                values = model_dict[variable]['layers'][layer]
                if not sortbylon:
                    arr=pd.DataFrame(np.vstack([hpar_list[hpar_idx],values]).T,columns =['lon', 'lat', 'pxw','values'])
                    arr = arr.sort_values(by=['lon','lat'])
                    values = arr['values'].values
                data_array[dep_indx[i],:,:] = np.reshape(values,
                                    (len(lat),len(lon)),order='F')
        else:
            data_array = xr.DataArray(np.zeros((len(lat),len(lon))),
                                      dims = ['latitude','longitude'],
                                      coords = [lat,lon])
            # if sorting is needed before reshaping
            values = model_dict[variable]['layers'][0]
            if not sortbylon:
                arr=pd.DataFrame(np.vstack([hpar_list[hpar_idx],values]).T,columns =['lon', 'lat', 'pxw','values'])
                arr = arr.sort_values(by=['lon','lat'])
                values = arr['values'].values
            data_array[:,:] = np.reshape(values,
                                    (len(lat),len(lon)),order='F')
        #-------------------------------------------------------------------------
        #add reference values at each depth as metadata to the Data_Array
        #-------------------------------------------------------------------------
        av_attrs = {}
        for keys in parser['parameters'][variable].keys():
            if (sys.version_info[:2] > (3, 0)):
                av_attrs[keys] = parser['parameters'][variable][keys]
            else:
                av_attrs[keys] = parser['parameters'][variable][keys].decode('utf-8')

        # read the 1D model if any of the reference values are not defined
        av_attrs['refmodel'] = parser['metadata']['refmodel']

        if len(data_array.shape) == 3: # if 3-D variable
            # get the variable values
            av_depth = deepcopy(data_array.depth.values)
            avgvalue = []; ifaverage =  True

            # calculate reference value if 1D model is read and absolute_unit is specified
            if ifread1D and 'absolute_unit' in av_attrs.keys():
                # get custom parameters that share a name with existing variables
                # e.g. vs_even6 for even degree variations up to degree 6
                ref1d.get_custom_parameter(variable)

                target_unit=av_attrs['absolute_unit']
                refvalue = ref1d.evaluate_at_depth(av_depth,parameter=variable).to(target_unit).magnitude
                av_attrs['refvalue'] = refvalue

            # loop over depths and get average values
            for _,depth in enumerate(av_depth):
                # select the appropriate map
                mapval = data_array.sel(depth=depth)
                # get the average, use an earlier evaluation of area if possible
                if ifaverage:
                    try:
                        globalav,area, _  = tools.meanxarray(mapval,area=area,pix_width = ds['pixel_width'])
                        avgvalue.append(globalav)
                    except:
                        print(traceback.format_exc())
                        warnings.warn('Could not read mean values for parameter '+variable)
                        ifaverage = False
            if ifaverage: av_attrs['average'] = np.array(avgvalue)
            av_attrs['start_depths'] = allstartdepths
            av_attrs['end_depths'] = allenddepths
        else:
            # get the average, use an earlier evaluation of area if possible
            try:
                globalav,area,_ = tools.meanxarray(data_array,area=area,pix_width = ds['pixel_width'])
                av_attrs['average'] = globalav
            except:
                warnings.warn('Could not read mean values for parameter '+variable)
            av_attrs['depth'] = float(av_attrs['depth'])

        #add Data_Array object to Data_Set
        data_array.attrs = av_attrs
        ds[variable] = data_array

    #Add overall attributes
    attrs = {}
    for key in parser['metadata'].keys():
        if (sys.version_info[:2] > (3, 0)):
            attrs[key] = parser['metadata'][key]
        else:
            attrs[key] = parser['metadata'][key].decode('utf-8')
    # add data variables as an attributes since pixel_width or other data may also be added
    attrs['parameters'] = variables
    ds.attrs = attrs

    # write to netcdf
    comp = {'zlib': True, 'complevel': complevel}
    encoding = {var: comp for var in ds.data_vars}
    if outfile != None:
        # change None to string since it cannot be stored in
        for key in ds.attrs.keys():
            if ds.attrs[key] is None: ds.attrs[key] = 'None'
        for var in ds.data_vars:
            for key in ds[var].attrs.keys():
                if ds[var].attrs[key] is None: ds[var].attrs[key] = 'None'
        ds.to_netcdf(outfile,engine=engine,encoding=encoding)

    return ds

def getLU2symmetric(insparse):
    """Get the full symmetric matrix from LU matrix

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    print(".... Converting from LU matrix to symmetric matrix")
    outsparse=insparse.tolil(copy=True)
    outsparse.setdiag(0.)
    outsparse=outsparse.tocsr()
    outsparse=outsparse+insparse.T
    return outsparse

def readResCov(infile: str, onlymetadata: bool = False):
    """Reads Resolution or Covariance matrix created by invwdata_pm64 with option -r.
    R=inv(ATA+DTD)ATA and the name of file is typically outmodel.Resolution.bin

    Parameters
    ----------
    infile : str
        Binary file containing Resolution or Covariance matrix
    onlymetadata : bool, optional
        Only read the metadata and ignore the elements of the matrix, by default False

    Returns
    -------
    refmdl : str
        Reference 1D model
    kerstr : str
        Kernel signifying the basis sets comprising radial and horizontal parameterization
    ntot : int
        Number of parameters or basis coefficients in the model
    indexrad1,indexrad2,indexhor1,indexhor2 : np.ndarray
        Arrays describing the radial and horizontal basis sets that each ATA element in `out` correponds to
    out : np.ndarray
        Elements of the sparse Resolution or Covariance matrix

    Raises
    ------
    IOError
        Input file has not been found in the directory
    ValueError
        Number of bytes in binary file does not match expected

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    #read all the bytes to indata
    if (not os.path.isfile(infile)): raise IOError("Filename (",infile,") does not exist")
    nbytes = os.path.getsize(infile)

    ii = 0 #initialize byte counter
    ifswp = '' # Assuem that byte order is not be swapped unless elat is absurdly high
    start_time = timeit.default_timer()

    with open(infile, "rb") as f:
        # preliminary metadata
        indata = f.read(4) # try to read iflag
        iflag = struct.unpack(ifswp+'i',indata)[0] # Read flag
        if iflag != 1:
            ifswp = '!' # swap endianness from now on
            iflag = struct.unpack(ifswp+'i',indata)[0]
            if iflag != 1: raise ValueError("Error: iflag != 1")
        refmdl = struct.unpack('80s',f.read(80))[0].strip().decode('utf-8')
        kerstr = struct.unpack('80s',f.read(80))[0].strip().decode('utf-8')
        ntot = struct.unpack(ifswp+'i',f.read(4))[0]
        ndtd = int(((ntot+1)*ntot)/2)

        # pre-allocate matrices
        indexrad1 = None if onlymetadata else np.zeros(ndtd,dtype=int)
        indexrad2 = None if onlymetadata else np.zeros(ndtd,dtype=int)
        indexhor1 = None if onlymetadata else np.zeros(ndtd,dtype=int)
        indexhor2 = None if onlymetadata else np.zeros(ndtd,dtype=int)
        out = None if onlymetadata else np.zeros(ndtd)

        if not onlymetadata:
            # Now start reading data
            for jj in range(ndtd):
                indexrad1[jj] = struct.unpack(ifswp+'i',f.read(4))[0]
                indexrad2[jj] = struct.unpack(ifswp+'i',f.read(4))[0]
                indexhor1[jj] = struct.unpack(ifswp+'i',f.read(4))[0]
                indexhor2[jj] = struct.unpack(ifswp+'i',f.read(4))[0]
                out[jj] = struct.unpack(ifswp+'d', f.read(8))[0]
    ii=168+ndtd*24
    if ii != nbytes: raise ValueError("Error: number of bytes read ",str(ii)," do not match expected ones ",str(nbytes))
    elapsed = timeit.default_timer() - start_time
    if not onlymetadata: print(".... read "+str(ndtd)+" rows for the Res or Cov matrix in "+str(round(elapsed/60*10)/10)+" min.")

    return refmdl, kerstr, ntot, indexrad1, indexrad2, indexhor1, indexhor2, out

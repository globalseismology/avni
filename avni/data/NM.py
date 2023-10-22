#!/usr/bin/env python

#######################################################################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

#####################  IMPORT STANDARD MODULES   ######################################

import os
import sys
import h5py
import struct
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from avni.models import Reference1D
from avni import constants,tools
import pint
import warnings

if sys.version_info[0] >= 3: unicode = str

####################       I/O ROUTINES     ######################################

def read_rts_catalog(infile, base_units = True):
    """reads a mode catalog in bimary rts format

    Input Parameters:
    ----------------

    base_units: convert from native units to base units in constants
    """
    #initialize variables
    metadata = {}
    modelvar = []
    units = []
    rad = []
    val_temp = {}
    ilev_flag = []

    # reference1D instance and the units class
    pint.PintType.ureg = constants.ureg
    ref1d = Reference1D()

    #open buffer
    f = open(infile, 'rb')
    nbytes = os.path.getsize(infile)
    cc = 0 #initialize byte counter

    #read header
    indata = f.read(4); cc += 4
    ifswp='!'
    iflag = struct.unpack(ifswp+'i',indata)[0]
    if iflag != 1:
        ifswp = ''
        iflag = struct.unpack(ifswp+'i',indata)[0]
        if iflag != 1: raise ValueError("iflag != 1. Something wrong with endianness")
    ref1d._description=struct.unpack('80s',f.read(80))[0].strip().decode('utf-8');cc += 80
    ref1d._name=struct.unpack('80s',f.read(80))[0].strip().decode('utf-8'); cc += 80
    ref1d._infile=struct.unpack('80s',f.read(80))[0].strip(); cc += 80
    ref1d.metadata['ifanis']=struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
    ref1d.metadata['ideck']=struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
    ref1d.metadata['ref_period']=struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
    ref1d._nlayers = struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
    ref1d.metadata['itopic'] = struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
    ref1d.metadata['itopoc'] = struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
    ref1d.metadata['itopmantle'] = struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
    ref1d.metadata['itopcrust'] = struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4

    #Get parameters and units
    for i in range(0,9):
        tmpstr = struct.unpack('20s',f.read(20))[0].strip().decode('utf-8'); cc += 20
        if len(tmpstr.split())==2: # try getting the first part of the string assuming second part is units
            modelvar.append(tmpstr.split()[0].lower())
            units.append(tmpstr.split()[1].lower())
            print(tmpstr.split()[1].lower())
        else:
            modelvar.append(tmpstr.lower())
            units.append('dimensionless')
    ref1d.metadata['parameters'] = modelvar
    ref1d.metadata['units'] = units

    for i in range(0,ref1d._nlayers):
        for var in modelvar:
            if i == 0: val_temp[var] = []
            val_temp[var].append(struct.unpack(ifswp+'f',f.read(4))[0]); cc += 4
        ilev_flag.append(struct.unpack(ifswp+'i',f.read(4))[0]); cc += 4

    #names=['rho','vpv','vsv','qkappa','qmu','vph','vsh','eta']
    #use units from the elas/anelas file
    names = tools.convert2nparray(ref1d.metadata['parameters'])
    units = tools.convert2nparray(ref1d.metadata['units'])
    fields=list(zip(names,units))

    # loop over names and call evaluate_at_depth
    # Create data array for converted to Panda array with units
    PA_ = pint.PintArray; temp_dict = {}
    for paraindx,param in enumerate(names):
        print(paraindx,param)

        temp_dict[param] = PA_(val_temp[param], dtype="pint["+units[paraindx]+"]")
        if '/' in param: # convert from fractions such as 1000/qmu to qmu
            frac = param.split('/')
            numerator = float(frac[0])
            param_new = frac[-1]
            # loop over names and check if there's /name; modify units if needed
            val_temp2 = np.divide(numerator,val_temp[param],out=np.zeros_like(val_temp[param]), where=val_temp!=0)
            temp_dict[param_new] = PA_(val_temp2, dtype="pint[1/"+units[paraindx]+"]")

    modelarr = pd.DataFrame(temp_dict)
    if base_units: # convert to base units
        for col in modelarr.columns: modelarr[col] = modelarr[col].pint.to_base_units()
    modelarr['depth'] = PA_((constants.R.magnitude - modelarr['radius'].pint.to(constants.R.units).data).tolist(), dtype = constants.R.units)

    # number of layers
    nlev_out=struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4

    # finish off
    ref1d._nlayers = len(modelarr['radius'])
    ref1d.data = modelarr
    ref1d._radius_max = max(ref1d.data['radius'])

    #scan modes
    nmod_par=struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
    mod_par_desc = []
    for i in range(0,nmod_par):
        mod_par_desc.append(struct.unpack('20s',f.read(20))[0].strip().decode('utf-8'))
        cc += 20
    nbytpmod=(nmod_par+6*nlev_out)*4
    nbytesmo=nbytpmod
    ibytfrst=cc
    if (nbytes-ibytfrst) % nbytpmod ==0:
        nmodes = int((nbytes-ibytfrst)/nbytpmod)
    else:
        #print('byte position before modes:'+str(ibytfrst))
        #print('bytes per mode:'+str(nbytpmod))
        raise ValueError('remaining number of bytes should be divible by number of bytes per mode')

    # revise the number of paramters for each mode
    while 'null' in mod_par_desc: mod_par_desc.remove('null')
    null_par = nmod_par - len(mod_par_desc)
    nmod_par = len(mod_par_desc)

    # start reading modes
    rts_catalog = {}
    nnmax=0
    llmax=0
    ommax=0.
    nmodrad=0
    nmodsph=0
    nmodtor=0
    #mode types RTS
    ST = ['S','T','S']

    for ii in range(0,nmodes):
        nn=struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
        itype = struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
        ll = struct.unpack(ifswp+'i',f.read(4))[0]; cc += 4
        omega = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        if ll==0 and itype==1: nmodrad=nmodrad+1
        if ll>0 and itype==2: nmodtor=nmodtor+1
        if ll>0 and itype==3: nmodsph=nmodsph+1
        nnmax=max(nnmax,nn)
        llmax=max(llmax,ll)
        ommax=max(ommax,omega)
        name=str(nn)+ST[itype-1]+str(ll)

        #print(name,nn,itype,ll,omega,(2.*np.pi)/omega)
        print(name,nn,itype,ll,omega,omega/(2.*np.pi))
        #if name=='0S2': pdb.set_trace()

        # read the mode properties based on putmodebuffer
        smallq = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        gvel = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        vacc = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        hacc = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        vdis = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        hdis = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        pot = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        pvel= 6371.*omega/(float(ll)+0.5)
        # get junk values for the rest since they are undefined with null
        junk = f.read(4*null_par); cc += 4*null_par

        #make rts_catalog dictionary for each mode
        rts_catalog[name] = {}
        rts_catalog[name]['nn'] = nn
        rts_catalog[name]['itype'] = itype
        rts_catalog[name]['ll'] = ll
        rts_catalog[name]['omega'] = omega
        rts_catalog[name]['gvel'] = gvel
        rts_catalog[name]['vacc'] = vacc
        rts_catalog[name]['hacc'] = vacc
        rts_catalog[name]['vdis'] = vdis
        rts_catalog[name]['hdis'] = hdis
        rts_catalog[name]['pot'] = pot
        rts_catalog[name]['pvel'] = pvel

        # eigenfunction info
        # figure this out later-basically these are values of U,V,P for various depths
        nbyts = nbytesmo-(nmod_par+null_par)*4
        junk = f.read(nbyts); cc += nbyts
        #buff=np.zeros(4); buff2=np.zeros(4)
        #for jj in range(4): buff[jj] = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        #for jj in range(4): buff2[jj] = struct.unpack(ifswp+'f',f.read(4))[0]; cc += 4
        #u,up,upp = splinenm(buff(1),buff(2),buff2(1),buff2(2),dr,hn,hn2,hn3)
        #v,vp,vpp = splinenm(buff(3),buff(4),buff2(3),buff2(4),dr,hn,hn2,hn3,v,vp,vpp)

    if cc != nbytes: raise ValueError("number of bytes read "+str(cc)+" do not match expected ones "+str(nbytes))

    return rts_catalog,ref1d

def write_modes_hdf(infile,table_name,absorption_band='None'):
    '''
    writes an hdf5 file based on the output of getgennomo(?)

    params:
    infile <str>: path to mode catalog produced by getgennomo
    table_name <str>: name of hdf5 modes table that will be written
    '''
    rts_catalog,ref1d = read_rts_catalog(infile)

    #open hdf5 file
    ds = h5py.File(table_name)
    ds.attrs['model'] = ref1d.name
    ds.attrs['ref_period'] = ref1d.metadata['ref_period']
    ds.attrs['absorption_band'] = absorption_band
    #TODO add any other model attributes?
    modeslib = ds.create_group('modes')

    #create groups for mode types
    ds.create_group('radial')
    ds.create_group('toroidal')
    ds.create_group('spheroidal')

    #setup lists for fundamental mode dispersion branches
    omega_radial = []
    pvel_radial = []
    gvel_radial = []
    omega_toroidal = []
    pvel_toroidal = []
    gvel_toroidal = []
    omega_spheroidal = []
    pvel_spheroidal = []
    gvel_spheroidal = []

    disp_curve_dict = {}
    disp_curve_dict['radial'] = {}
    disp_curve_dict['toroidal'] = {}
    disp_curve_dict['spheroidal'] = {}

    for mode_name in rts_catalog:

        print('writing mode: {} to mode table'.format(mode_name))

        itype = rts_catalog[mode_name]['itype'] #mode type
        ll = rts_catalog[mode_name]['ll'] #angular order
        nn = rts_catalog[mode_name]['nn'] #radial order

        #create a new group for each radial order
        if itype == 1:
            try:
                ds['radial'].create_group(str(nn))
                disp_curve_dict['radial'][str(nn)] = {}
                disp_curve_dict['radial'][str(nn)]['omega'] = []
                disp_curve_dict['radial'][str(nn)]['gvel'] = []
                disp_curve_dict['radial'][str(nn)]['pvel'] = []
            except:
                warnings.warn('Could not create group for radial modes')
                pass

            ds['radial'][str(nn)].create_group(mode_name)
            ds['radial'][str(nn)][mode_name].create_dataset('omega',data=rts_catalog[mode_name]['omega'])
            ds['radial'][str(nn)][mode_name].create_dataset('pvel',data=rts_catalog[mode_name]['pvel'])
            ds['radial'][str(nn)][mode_name].create_dataset('gvel',data=rts_catalog[mode_name]['gvel'])
            ds['radial'][str(nn)][mode_name].create_dataset('vacc',data=rts_catalog[mode_name]['vacc'])
            ds['radial'][str(nn)][mode_name].create_dataset('hacc',data=rts_catalog[mode_name]['hacc'])
            ds['radial'][str(nn)][mode_name].create_dataset('vdis',data=rts_catalog[mode_name]['vdis'])
            ds['radial'][str(nn)][mode_name].create_dataset('hdis',data=rts_catalog[mode_name]['hdis'])
            ds['radial'][str(nn)][mode_name].create_dataset('pot',data=rts_catalog[mode_name]['pot'])

            disp_curve_dict['radial'][str(nn)]['omega'].append(rts_catalog[mode_name]['omega'])
            disp_curve_dict['radial'][str(nn)]['pvel'].append(rts_catalog[mode_name]['pvel'])
            disp_curve_dict['radial'][str(nn)]['gvel'].append(rts_catalog[mode_name]['gvel'])

        elif itype == 2:
            try:
                ds['toroidal'].create_group(str(nn))
                disp_curve_dict['toroidal'][str(nn)] = {}
                disp_curve_dict['toroidal'][str(nn)]['omega'] = []
                disp_curve_dict['toroidal'][str(nn)]['gvel'] = []
                disp_curve_dict['toroidal'][str(nn)]['pvel'] = []
            except:
                warnings.warn('Could not create group for toroidal modes')
                pass

            ds['toroidal'][str(nn)].create_group(mode_name)
            ds['toroidal'][str(nn)][mode_name].create_dataset('omega',data=rts_catalog[mode_name]['omega'])
            ds['toroidal'][str(nn)][mode_name].create_dataset('pvel',data=rts_catalog[mode_name]['pvel'])
            ds['toroidal'][str(nn)][mode_name].create_dataset('gvel',data=rts_catalog[mode_name]['gvel'])
            ds['toroidal'][str(nn)][mode_name].create_dataset('vacc',data=rts_catalog[mode_name]['vacc'])
            ds['toroidal'][str(nn)][mode_name].create_dataset('hacc',data=rts_catalog[mode_name]['hacc'])
            ds['toroidal'][str(nn)][mode_name].create_dataset('vdis',data=rts_catalog[mode_name]['vdis'])
            ds['toroidal'][str(nn)][mode_name].create_dataset('hdis',data=rts_catalog[mode_name]['hdis'])
            ds['toroidal'][str(nn)][mode_name].create_dataset('pot',data=rts_catalog[mode_name]['pot'])

            disp_curve_dict['toroidal'][str(nn)]['omega'].append(rts_catalog[mode_name]['omega'])
            disp_curve_dict['toroidal'][str(nn)]['pvel'].append(rts_catalog[mode_name]['pvel'])
            disp_curve_dict['toroidal'][str(nn)]['gvel'].append(rts_catalog[mode_name]['gvel'])

        elif itype == 3:
            try:
                ds['spheroidal'].create_group(str(nn))
                disp_curve_dict['spheroidal'][str(nn)] = {}
                disp_curve_dict['spheroidal'][str(nn)]['omega'] = []
                disp_curve_dict['spheroidal'][str(nn)]['gvel'] = []
                disp_curve_dict['spheroidal'][str(nn)]['pvel'] = []
            except:
                warnings.warn('Could not create group for spheroidal modes')
                pass

            ds['spheroidal'][str(nn)].create_group(mode_name)
            ds['spheroidal'][str(nn)][mode_name].create_dataset('omega',data=rts_catalog[mode_name]['omega'])
            ds['spheroidal'][str(nn)][mode_name].create_dataset('pvel',data=rts_catalog[mode_name]['pvel'])
            ds['spheroidal'][str(nn)][mode_name].create_dataset('gvel',data=rts_catalog[mode_name]['gvel'])
            ds['spheroidal'][str(nn)][mode_name].create_dataset('vacc',data=rts_catalog[mode_name]['vacc'])
            ds['spheroidal'][str(nn)][mode_name].create_dataset('hacc',data=rts_catalog[mode_name]['hacc'])
            ds['spheroidal'][str(nn)][mode_name].create_dataset('vdis',data=rts_catalog[mode_name]['vdis'])
            ds['spheroidal'][str(nn)][mode_name].create_dataset('hdis',data=rts_catalog[mode_name]['hdis'])
            ds['spheroidal'][str(nn)][mode_name].create_dataset('pot',data=rts_catalog[mode_name]['pot'])

            disp_curve_dict['spheroidal'][str(nn)]['omega'].append(rts_catalog[mode_name]['omega'])
            disp_curve_dict['spheroidal'][str(nn)]['pvel'].append(rts_catalog[mode_name]['pvel'])
            disp_curve_dict['spheroidal'][str(nn)]['gvel'].append(rts_catalog[mode_name]['gvel'])
        else:
            raise ValueError('itype = {} not recognized'.format(itype))

    #add dispersion curves to attributes of a mode overtone
    for item in disp_curve_dict['toroidal'].keys():
        ds['toroidal'][item].attrs['omega'] = np.array(disp_curve_dict['toroidal'][item]['omega'])
        ds['toroidal'][item].attrs['pvel'] = np.array(disp_curve_dict['toroidal'][item]['pvel'])
        ds['toroidal'][item].attrs['gvel'] = np.array(disp_curve_dict['toroidal'][item]['gvel'])
    for item in disp_curve_dict['spheroidal'].keys():
        ds['spheroidal'][item].attrs['omega'] = np.array(disp_curve_dict['spheroidal'][item]['omega'])
        ds['spheroidal'][item].attrs['pvel'] = np.array(disp_curve_dict['spheroidal'][item]['pvel'])
        ds['spheroidal'][item].attrs['gvel'] = np.array(disp_curve_dict['spheroidal'][item]['gvel'])

    #TODO write model attributes and max degree / overtone for S and T modes

def get_mode_attribute(table,mode_type,radial_order,angular_order,attribute):
    '''
    returns the eigenfrequency of a single mode

    params:
    table <str>: path to hdf5 format modes table
    mode_type <str>: either "spheroidal", "toroidal", or "radial"
    radial_order <int>: radial order (ie. overtone number)
    angular_order <int>: angular order
    attribute: characteristics of a mode such as pvel,gvel,omega

    returns:
    value: value of attribute queried
    '''

    table = h5py.File(table,'r')

    if mode_type=='S':
        mode_type='spheroidal'
    elif mode_type=='T':
        mode_type='toroidal'
    elif mode_type=='R':
        mode_type='radial'

    if mode_type == 'radial':
        mode_name = '{}R{}'.format(radial_order,angular_order)
    elif mode_type == 'toroidal':
        mode_name = '{}T{}'.format(radial_order,angular_order)
    elif mode_type == 'spheroidal':
        mode_name = '{}S{}'.format(radial_order,angular_order)

    value = table[mode_type][str(radial_order)][mode_name][attribute][()]

    return value


def get_mode_freq(table,mode_type,radial_order,angular_order,freq_units='mhz'):
    '''
    returns the eigenfrequency of a single mode

    params:
    table <str>: path to hdf5 format modes table
    mode_type <str>: either "spheroidal", "toroidal", or "radial"
    radial_order <int>: radial order (ie. overtone number)
    angular_order <int>: angular order

    returns:
    freq: eigen frequency of the mode in freq_units
    '''

    omega = get_mode_attribute(table,mode_type,radial_order,angular_order,'omega')

    if freq_units.lower() == 'hz':
        freq = omega/(2.*np.pi)
    elif freq_units.lower() == 'mhz':
        freq = (omega/(2.*np.pi)) * 1000.
    elif freq_units.lower() == 'rad/s':
        freq = omega

    return freq

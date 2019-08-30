#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int,list,tuple

#importing standard modules
import os
import glob
import fortranformat as ff #reading/writing fortran formatted text
import numpy as np
import xarray as xr
import pdb

############             input REM          modules ######
from .common import get_fullpath
from .bases import eval_ylm
###########################  ANALYSIS SUBROUTINES ####################

def getdepthsfolder(folder='.',extension='.epix',delimiter='.'):
    """"Get list of depths from filenames to iterate through"""
    depths = []
    folder = get_fullpath(folder)
    onlyfiles = glob.glob(folder+ '/*'+extension)
    for name in onlyfiles:
        name = name.split(folder)[1]
        depths.append(int(name[name.index(delimiter)+1:name.rindex(extension)]))
    depths.sort()
    return depths

def rdswpsh(filename):
    """Code to get spherical harmonic coefficients in the ylm normalization. shmatrix is the linear array
    with the ylm normalization like in PM's codes. """

    if (not os.path.isfile(filename)): raise IOError("Filename ("+filename+") does not exist")

    lineformat = ff.FortranRecordReader("(i4,i4,2e13.5)")
    comments=[]
    shmatrix=[]
    metadata={}
    for line in open(filename):
        if line.startswith("#"):
            if ':' in line:
                field = line.split(':')[0].split('#')[1].strip()
                metadata[field] = line.split(':')[1].split('\n')[0].strip()
            else:
                comments.append(line.split('\n')[0].strip())
        else:
            tempshrow = lineformat.read(line)
            # Change based on ylm normalization
            l=tempshrow[0];m=tempshrow[1]
            if m != 0:
                tempshrow[2] =  2. * tempshrow[2]
                tempshrow[3] = -2. * tempshrow[3]
            shmatrix.append(tempshrow)

    # Convert to numpy array
    namelist = ['l','m','cos','sin']
    formatlist = ['i4','i4','f8','f8']
    dtype = dict(names = namelist, formats=formatlist)
    shmatrix = np.array([tuple(x) for x in shmatrix],dtype=dtype)
    return shmatrix, metadata, comments


def wrswpsh(filename,shmatrix,metadata=None,comments=None, lmax=None):
    """Code to write spherical harmonic coefficients from the ylm normalization. shmatrix is the linear array
    with the ylm normalization like in PM's codes.

    Parameters
    ----------

    filename : Name of the file containing four columns
              (degree, order, cosin, sin)

    metadata: metadata fields from input fields if specified.
              default : {'FORMAT':'0'}

    comments: all other comments except lines containing metadata
    """
    # default value
    if metadata is None: metadata = {'FORMAT':'0'}

    #combine headers
    printstr=[]
    for key in sorted(metadata.keys()): printstr.append('#'+key+':'+metadata[key]+'\n')
    if comments is not None:
        for comment in comments: printstr.append(comment+'\n')
    header_linem0 = ff.FortranRecordWriter('(i4,i4,e13.5)')
    header_line = ff.FortranRecordWriter('(i4,i4,2e13.5)')
    for ii,degree in enumerate(shmatrix['l']):
        order = shmatrix['m'][ii]
        if lmax == None: lmax = degree + 1000 # include all degrees if None
        if degree <= lmax:
            if order == 0:
                arow=header_linem0.write([degree,order,shmatrix['cos'][ii]])
            else:
                # convert to ylm normalization on the fly
                arow=header_line.write([degree,order,0.5*shmatrix['cos'][ii],-0.5*shmatrix['sin'][ii]])
            printstr.append(arow+'\n')

    # write out the file
    f=open(filename,'w')
    f.writelines(printstr)
    f.close()
    print("..... written "+filename)
    return

def get_coefficients(shmatrix,lmin=None,lmax=None):
    if lmin == None: lmin = min(shmatrix['l'])
    if lmax == None: lmax = max(shmatrix['l'])
    if lmin == 0:
        coeff = np.zeros((lmax+1)**2)
    else:
        coeff = np.zeros((lmax+1)**2-lmin**2)
    iloop = 0
    for degree in np.arange(lmin,lmax+1):
        select = (shmatrix['l']==degree) & (shmatrix['m']==0)
        coeff[iloop] = shmatrix['cos'][select]
        iloop += 1
        for order in np.arange(1,degree+1):
            select = (shmatrix['l']==degree) & (shmatrix['m']==order)
            coeff[iloop] = shmatrix['cos'][select]
            iloop += 1
            coeff[iloop] = shmatrix['sin'][select]
            iloop += 1
    return coeff

def swp_to_xarray(shmatrix,grid=10,lmax=None):
    latitude = np.arange(-90+grid/2., 90,grid)
    longitude = np.arange(0+grid/2., 360,grid)
    if lmax == None: lmax = max(shmatrix['l'])
    coeff = get_coefficients(shmatrix,lmax=lmax)
    values = eval_ylm(latitude,longitude,lmaxhor=lmax,grid=True,norm='ylm',weights=coeff)
    nlat = len(latitude)
    nlon = len(longitude)

    # index column by longitude
    values = values.reshape((nlat,nlon),order='C')

    # get the grid sizes stored
    outarr = xr.DataArray(values, dims = ['latitude','longitude'],
                        coords = [latitude,longitude])
    return outarr


def calcshpar2(shmatrix):
    """This calculates the roughness and rms of a model file from shmatrix read using rdswpsh"""

    # convert from shmatrix to shmod
    lmod = max(shmatrix['l'])
    icount = 0
    shmod = []
    for ii in np.arange(len(shmatrix['l'])):
        if shmatrix['m'][ii] == 0:
            shmod.append(shmatrix['cos'][ii])
            icount = icount + 1
        else:
            shmod.append(shmatrix['cos'][ii])
            shmod.append(shmatrix['sin'][ii])
            icount = icount + 2
    if icount != (lmod+1)**2: raise ValueError("icount != (lmod+1)**2")

#     loop over all degrees, fixing the normalization from ylm on the
#     fly
    i=0
    powerarr = np.zeros(lmod+1)
    for l in np.arange(0,lmod+1):
        power=0.
        for m in np.arange(0,l+1):
            if m == 0:
                coeff=shmod[i]/np.sqrt(4.*np.pi)
                power=power+coeff**2
                i=i+1
            else:
                coeff=np.sqrt(2.)*0.5*shmod[i]/np.sqrt(4.*np.pi)
                power=power+coeff**2
                coeff=-np.sqrt(2.)*0.5*shmod[i+1]/np.sqrt(4.*np.pi)
                power=power+coeff**2
                i=i+2
        powerarr[l]=power
    powerall=0.
    powerall0=0.
    gradsquared=0.
    xlaplacesquared=0.
#
    average=shmod[0]/np.sqrt(4.*np.pi)
    for l in np.arange(0,lmod+1):
        yy=powerarr[l]
        if l != 0:
            powerall=powerall+yy
            gradsquared=gradsquared+yy*float(l)*float(l+1)
            xlaplacesquared=xlaplacesquared+yy*float(l)*float(l+1)*float(l)*float(l+1)
        powerall0=powerall0+yy
    powerall=np.sqrt(powerall)
    powerall0=np.sqrt(powerall0)
#
# powerall0 is the rms of all ls, including l=0
# powerall is the rms of all l>0
    rms=powerall
    roughness=np.sqrt(gradsquared)/rms
#    roughness2=sqrt(gradsquared)/rms
#    roughness1=0.

    return average,rms,roughness,powerarr

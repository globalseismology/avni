#!/usr/bin/env python

# python 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import *

#importing standard modules
import os
import glob
import pdb
import fortranformat as ff #reading/writing fortran formatted text
import numpy as np
############             input REM          modules ######
from .common import get_fullpath
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


def rdswpsh(filename,maxl=360):
    """Code to get spherical harmonic coefficients in the ylm normalization. shmatrix is the linear array 
    with the ylm normalization like in PM's codes. """
    
    if (not os.path.isfile(filename)): raise IOError("Filename ("+filename+") does not exist")
        
    lineformat = ff.FortranRecordReader("(i4,i4,2e13.5)")
    comments=[]
    shmatrix=[]
    for line in open(filename):
        if line.startswith("#"): 
            comments.append(line)
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
    return shmatrix, comments

def wrswpsh(filename,shmatrix,comments, maxl=360):
    """Code to write spherical harmonic coefficients from the ylm normalization. shmatrix is the linear array 
    with the ylm normalization like in PM's codes. """
    
    header_linem0 = ff.FortranRecordWriter('(i4,i4,e13.5)')        
    header_line = ff.FortranRecordWriter('(i4,i4,2e13.5)')
    printstr = list(comments)    
    for ii in np.arange(len(shmatrix['l'])):
        if shmatrix['l'][ii] <= maxl:
            if shmatrix['m'][ii] == 0:
                arow=header_linem0.write([shmatrix['l'][ii],shmatrix['m'][ii],shmatrix['cos'][ii]])
            else:
                # convert to ylm normalization on the fly
                arow=header_line.write([shmatrix['l'][ii],shmatrix['m'][ii],0.5*shmatrix['cos'][ii],-0.5*shmatrix['sin'][ii]])
            
            printstr.append(arow+'\n')
    
    f=open(filename,'w')
    f.writelines(printstr)
    f.close()
    print("..... written "+filename)
    return
 
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

#!/usr/bin/env python
"""This script/module contains routines that are used to analyze/visualize the data sets 
in the standard REM3D format. 
Author: Raj Moulik, 2018"""

#####################  IMPORT STANDARD MODULES   ######################################   

import sys,os
import argparse #parsing arguments
import glob # pattern matching for finding files
import numpy as np #for numerical analysis
from scipy import stats #for stats analysis
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output and compare times
import pdb    #for the debugger pdb.set_trace()
from math import pi
import matplotlib.pyplot as plt
import fortranformat as ff #reading/writing fortran formatted text

#######################################################################################


def readREM3DSWformat(file):
    """Reads the REM3D format for analysis and plotting"""

    if (not os.path.isfile(file)): sys.exit("Filename ("+file+") does not exist")
    namelist=['overtone','peri','typeiorb','cmtname','eplat','eplon','cmtdep', \
              'stat','stlat','stlon','distkm','refphase','delobsphase', \
              'delerrphase','delpredphase']
    formatlist=['i','f8','a2','a15','f8','f8',  \
                'f8','a8','f8','f8','f8','f8',  \
                'f8','f8','f8']

    dtype = dict(names = namelist, formats=formatlist)

    # checks for CITE and SHORTCITE comments
    comments=[]
    for line in open(file):
        if line.startswith("#"): comments.append(line)
    reference='None'
    for line in open(file):
        if line.startswith("#CITE:"): 
            reference=line[6:].rstrip('\n')
            break
    if reference=='None': print "Reference (e.g. #CITE: Ekstrom, 2011) should be defined as a comment in "+file ; sys.exit(2)
    shortref='None'
    for line in open(file):
        if line.startswith("#SHORTCITE:"): 
            shortref=line[11:].rstrip('\n')
            break
    if shortref=='None': print "Short reference name (e.g. #SHORTCITE: GDM52) should be defined as a comment in "+file ; sys.exit(2)

    SWdata = np.genfromtxt(file, dtype=dtype,names=namelist,comments="#")
    
    return SWdata,comments,reference,shortref


def writeREM3DSWformat(filename,SWdata,comments):
    """Writes the REM3D format for analysis and plotting"""
    
    header_line = ff.FortranRecordWriter('(i2,1x,f6.1,1x,a2,1x,a15,1x,f8.3,1x,f9.3,1x, \
                                           f6.1,1x,a8,1x,f8.3,1x,f9.3,1x,f9.1,1x,f10.3,1x, \
                                           f10.3,1x,f10.3,1x,f10.3)')
    printstr = list(comments)
    for ii in np.arange(len(SWdata['overtone'])):
        arow=header_line.write([SWdata['overtone'][ii],SWdata['peri'][ii],SWdata['typeiorb'][ii],SWdata['cmtname'][ii].ljust(15),SWdata['eplat'][ii],SWdata['eplon'][ii],SWdata['cmtdep'][ii],SWdata['stat'][ii].ljust(8),SWdata['stlat'][ii],SWdata['stlon'][ii],SWdata['distkm'][ii],SWdata['refphase'][ii],SWdata['delobsphase'][ii],SWdata['delerrphase'][ii],SWdata['delpredphase'][ii]])
        printstr.append(arow+'\n')
    f=open(filename,'w')
    f.writelines(printstr)
    f.close()
    print "....written "+str(len(SWdata))+" observations to "+filename
    return

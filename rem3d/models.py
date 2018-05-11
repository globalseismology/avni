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
import fortranformat as ff #reading/writing fortran formatted text

#####################

def readepixfile(filename):
	"""Read .epix file format."""

	currentdir=os.getcwd()
	try: 
		f = open(filename, 'r')
	except IOError:
		print "File ("+filename+") does not exist in the current directory - "+currentdir
		sys.exit(2)
			
	epixarr=np.genfromtxt(filename, dtype=None,comments="#",names=['lat','lon','pixsize','val'])

	return epixarr
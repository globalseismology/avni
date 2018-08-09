#!/usr/bin/env python

from __future__ import absolute_import
import re
import os
import sys
import subprocess
import glob

# import numpy
#---------------------------------------------------------------
try:
    import numpy
except:
    raise ImportError( 'REM3D requires Numpy 1.7 or later.' )

npv = numpy.__version__.split('.')
if int(npv[0]) != 1:
    raise ImportError( 'REM3D requires Numpy 1.7 or later.' )
if int(npv[1]) < 7:
    raise ImportError( 'REM3D requires Numpy 1.7 or later.' )

# write short description
#--------------------------------------------------------------------------
description = 'a modeling and analysis toolkit for reference Earth ' + \
    'datasets and tomographic models.'


# puts the contents of the README file in the variable long_description
#--------------------------------------------------------------------------
try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md').read()

# call setup
#--------------------------------------------------------------------------

versionstuff = dict(
    re.findall("(.+) = '(.+)'\n", open('rem3d/version.py').read()))


# Tried to use setuptools in order to check dependencies.
# if the system does not have setuptools, fall back on
# distutils.
# Can only use numpy distutils with distutils as fortran codes 
# are not compiles otherwise
# Build the f2py fortran extension
# --------------------------------
from os.path import join
from numpy.distutils.core import Extension
from numpy.distutils.core import setup
 
metadata = dict(name = 'rem3d',
                version=versionstuff['version'],
                description=description,
                long_description = long_description,
                long_description_content_type='text/markdown',
                url='http://www.rem3d.org',
                author = 'Pritwiraj Moulik',
                author_email='pritwiraj.moulik@gmail.com',
                license='GPL',
                packages = ['rem3d'],
                install_requires=['pandas','fortranformat','joblib','progressbar',
                'netCDF4'],
                data_files=[('rem3d', ['README.md'])],
                keywords = ['earth-science','earth-observation','earthquake',
                'earth','earthquake-data','geology','geophysics',
                'geophysical-inversions','seismology','seismic-inversion',
                'seismic-waves','seismic-tomography','mineral',
                'geochemistry','geodesy','physics','modeling','modeling-tool',
                'model','geodynamics'],
                classifiers=[
                'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                'Programming Language :: Fortran',
                'Programming Language :: Python',
                'Intended Audience :: Science/Research',
                'Topic :: Education',
                'Natural Language :: English',
                ],
                )

setup(**metadata)

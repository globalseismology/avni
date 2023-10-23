#!/usr/bin/env python

from __future__ import absolute_import
import re
import os
#import sys
#import subprocess
#import glob

# import numpy
#---------------------------------------------------------------
try:
    import numpy
except:
    raise ImportError( 'AVNI requires Numpy 1.7 or later.' )

npv = numpy.__version__.split('.')
if int(npv[0]) != 1:
    raise ImportError( 'AVNI requires Numpy 1.7 or later.' )
if int(npv[1]) < 7:
    raise ImportError( 'AVNI requires Numpy 1.7 or later.' )

import numpy.distutils.fcompiler

# get F90 environment variable
#----------------------------------------
# F90 = os.getenv("F90")
F90 = "gfortran"

# raise error if F90 not defined
# !! comment out this if statement for manual install !!
#------------------------------------------------------------
if F90 == None or F90 == "":
    l1 = 'AVNI requires environment variable F90 to be set. \n '
    l2 = 'Please set to one of {"ifort", "gfortran"}'
    raise RuntimeError( l1 + l2 )


# specialize for different compilers
#------------------------------------------------------------
if F90 == "ifort":
    f90_flags = ["-openmp", "-fPIC", "-xHost", "-O3", "-ipo",
                 "-funroll-loops", "-heap-arrays", "-mcmodel=medium"]
    omp_lib = ["-liomp5"]

elif F90 == "gfortran":
    f90_flags = ["-libmil","-fopenmp", "-fPIC", "-O3", "-fbounds-check",
                 "-ffixed-line-length-none","-std=legacy"]
    omp_lib = ["-lgomp"]

elif F90 == "f90":
    f90_flags = ["-fopenmp", "-fPIC", "-O3", "-library=sunperf","-xopenmp"]
    omp_lib = [""]

elif F90 in ["pgfortran", "pgf90", "pgf95"]:
    f90_flags = ["-mp"]
    omp_lib = [""]

else:
    l1 = "F90 = " + F90 + ". \n"
    l2 = "Environment variable F90 not recognized.  \n"
    raise RuntimeError( l1 + l2 )

# use old version of memcpy
# https://snorfalorpagus.net/blog/2016/07/17/compiling-python-extensions-for-old-glibc-versions/
# os.environ['CFLAGS']="-I. -include docs/.glibc_version_fix.h"

# for manual install comment out the above section and define
# the variables f90_flags and omp_lib below
#------------------------------------------------------------
#f90_flags =
#omp_lib =


# discover fortran compiler
#---------------------------------------------------------------
#cmnd = 'f2py -c --help-fcompiler'
#fcstr = subprocess.check_output( ["f2py", "-c", "--help-fcompiler"] )
#ii = fcstr.index('Fortran compilers found:')
#ff = fcstr.index('Compilers available for this platform, but not found:')
#fcstr_cut = fcstr[ii:ff].split('\n')
#fcs = [ fcstr_cut[i].strip() for i in range(1,len(fcstr_cut)-1) ]
#print fcstr_cut
#print len(fcstr_cut)
#print fcs

# setup fortran 90 extension
#---------------------------------------------------------------------------

f90_dir='avni/f2py'
packagelist=['avni','avni.data','avni.models','avni.tools',
             'avni.mapping','avni.plots']
for module in os.listdir(f90_dir): packagelist.append('avni.f2py.'+module)


# write short description
#--------------------------------------------------------------------------
description = 'Analysis and Visualization toolkit for plaNetary Inferences ' + \
    'provides web-based and backend code access to tools, techniques, models and data.'


# puts the contents of the README file in the variable long_description
#--------------------------------------------------------------------------
file = open('README.md', 'r')
Lines = file.readlines()
long_description = '\n\n '
for string in Lines:
    if '<img' not in string: long_description += string

# call setup
#--------------------------------------------------------------------------

versionstuff = dict(
    re.findall("(.+) = '(.+)'\n", open('avni/version.py').read()))


# Tried to use setuptools in order to check dependencies.
# if the system does not have setuptools, fall back on
# distutils.
# Can only use numpy distutils with distutils as fortran codes
# are not compiles otherwise
#
# Build the f2py fortran extension
# --------------------------------
from os.path import join
import setuptools # added to make the Markdown long_description work for PyPi
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

# Use this if you need import avni.module for every module folder
# extf = [Extension(name='avni.'+module,
#                 sources = [join(f90_dir,module,f) for f in os.listdir(join(f90_dir,module)) if f.endswith('.f')],
#                 extra_f77_compile_args = f90_flags,
#                 extra_f90_compile_args = f90_flags,
#                 extra_link_args = omp_lib)
#             for module in os.listdir(f90_dir)]
#
# Use this if you need a single module for all subroutines import avni.f2py
sourcefiles = []
for path,_,filelist in os.walk(join(f90_dir)):
    for f in filelist:
        if f.endswith('.f'): sourcefiles.append(join(path,f))
extf = [Extension(name='avni.f2py',
                sources = sourcefiles,
                extra_f77_compile_args = f90_flags,
                extra_f90_compile_args = f90_flags,
                extra_link_args = omp_lib)]

def parse_requirements_file(fname):
    requirements = list()
    with open(fname, 'r') as fid:
        for line in fid:
            req = line.strip()
            if req.startswith('#'):
                continue
            # strip end-of-line comments
            req = req.split('#', maxsplit=1)[0].strip()
            requirements.append(req)
    return requirements

# data_dependencies is empty, but let's leave them so that we don't break
# people's workflows who did `pip install avni[all]`
install_requires = parse_requirements_file('requirements_base.txt')
all_requires = parse_requirements_file('requirements_extra_pip.txt')

metadata = dict(name = 'avni',
                version=versionstuff['version'],
                description=description,
                long_description_content_type='text/markdown',
                long_description = long_description,
                url='https://avni.globalseismology.org',
                author = 'Pritwiraj Moulik',
                author_email='moulik@caa.columbia.edu',
                license='GPL',
                packages = packagelist,
                ext_modules = extf,
                # Notes on why numpy and setuptools version are needed
                # https://numpy.org/devdocs/reference/distutils_status_migration.html
                data_files=[('avni', ['README.md','CODE_OF_CONDUCT.md','CONTRIBUTING.md',
                                      'requirements_base.txt','requirements.txt',
                                      'requirements_extra_pip.txt']),
                ('avni/config',['avni/config/attributes.ini',
                'avni/config/planets.ini','avni/config/units.ini',
                'avni/config/f2py.ini'])],
                keywords = ['earth-science','earth-observation','earthquake',
                'earth','earthquake-data','geology','geophysics',
                'geophysical-inversions','seismology','seismic-inversion',
                'seismic-waves','seismic-tomography','mineral',
                'geochemistry','geodesy','physics','modeling','modeling-tool',
                'model','geodynamics'],
                classifiers=[
                'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                'Programming Language :: Fortran',
                'Programming Language :: Python :: 3',
                'Intended Audience :: Science/Research',
                'Topic :: Education',
                'Natural Language :: English',
                ],
                project_urls={
                    'Documentation': 'https://avni.globalseismology.org',
                    'Source': 'https://github.com/globalseismology/avni',
                    'Tracker': 'https://github.com/globalseismology/avni/issues',
                },
                platforms='any',
                python_requires='>=3.7',
                install_requires=install_requires,
                extras_require={
                    'all': all_requires,
                },
                )

setup(**metadata)

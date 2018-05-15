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

import numpy.distutils.fcompiler

# get F90 environment variable
#----------------------------------------
F90 = os.getenv("F90")


# raise error if F90 not defined
# !! comment out this if statement for manual install !!
#------------------------------------------------------------
if F90 == None or F90 == "":
    l1 = 'REM3D requires environment variable F90 to be set. \n '
    l2 = 'Please set to one of {"ifort", "gfortran"}'
    raise RuntimeError( l1 + l2 )


# specialize for different compilers
#------------------------------------------------------------
if F90 == "ifort":
    f90_flags = ["-openmp", "-fPIC", "-xHost", "-O3", "-ipo",
                 "-funroll-loops", "-heap-arrays", "-mcmodel=medium"]
    omp_lib = ["-liomp5"]

elif F90 == "gfortran":
    f90_flags = ["-fopenmp", "-fPIC", "-O3", "-fbounds-check",
                 "-mtune=native"]
    omp_lib = ["-lgomp"]

else:
    l1 = "F90 = " + F90 + ". \n"
    l2 = "Environment variable F90 not recognized.  \n"
    raise RuntimeError( l1 + l2 )


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

f90_dir='rem3d/f2py'
packagelist=['rem3d']
for module in os.listdir(f90_dir): packagelist.append('rem3d.f2py.'+module)


# write short description
#--------------------------------------------------------------------------
description = 'a modeling and analysis toolkit for reference Earth ' + \
    'datasets and tomographic models.'


# puts the contents of the README file in the variable long_description
#--------------------------------------------------------------------------
with open('README.md') as file:
    long_description = '\n\n ' + file.read()


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
 
extf = [Extension(name='rem3d.'+module,
                sources = [join(f90_dir,module,f) for f in os.listdir(join(f90_dir,module)) if f.endswith('.f')],
                extra_f90_compile_args = f90_flags,
                extra_link_args = omp_lib)
            for module in os.listdir(f90_dir)]


metadata = dict(name = 'rem3d',
                version=versionstuff['version'],
                description=description,
                long_description = long_description,
                url='http://www.rem3d.org',
                author = 'Pritwiraj Moulik',
                author_email='pritwiraj.moulik@gmail.com',
                license='GPL',
                packages = packagelist,
                ext_modules = extf,
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

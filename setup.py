from __future__ import absolute_import
import re
import os
import sys
import subprocess


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

import numpy.distutils.core
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

elif F90 in ["pgfortran", "pgf90", "pgf95"]:
    f90_flags = ["-mp"]
    omp_lib = [""]

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

f90_fnames = [
    'types.f90',
    'utils.f90',
    'zhang_jin.f90',
    'zhang_jin_f.f90',
    'slatec.f90',
    'special_functions.f90',
    'legendre_polynomial.f90',
    'm_mrgrnk.f90',
    'physical_constants.f90',
    'geometry.f90',
    'hui_gnedin_97.f90',
    'hm12.f90',
    'chem_cool_rates.f90',
    'ion_solver.f90',
    'verner_96.f90',
    'photo_xsections.f90',
    'spectra.f90',
    'source_point.f90',
    'source_plane.f90',
    'source_background.f90',
    'sphere_base.f90',
    'sphere_stromgren.f90',
    'sphere_bgnd.f90',
    'slab_base.f90',
    'slab_plane.f90',
    'slab_bgnd.f90',
    ]



f90_paths = []
for fname in f90_fnames:
    f90_paths.append( 'rem3d/f2py/' + fname )




extf = numpy.distutils.core.Extension(
    name = 'rem3d-f90',
    sources = f90_paths,
    extra_f90_compile_args = f90_flags,
    extra_link_args = omp_lib,
    )


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

metadata = dict(name = 'rem3d',
                version=versionstuff['version'],
                description=description,
    			long_description = long_description,
                url='http://www.rem3d.org',
    			author = 'Pritwiraj Moulik',
                author_email='pritwiraj.moulik@gmail.com',
                license='GPL',
                long_description='REM3D is a Python library for reference Earth datasets and tomographic models.',
    			packages = ['rem3d'],
                package_data={'rem3d': ['*.f90','data/input_*/*']},
                ext_modules = [extf],
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

# Try to use setuptools in order to check dependencies.
# if the system does not have setuptools, fall back on
# distutils.
try:
    from setuptools import setup
    metadata['install_requires'] = ['numpy', 'matplotlib', 'scipy']
except ImportError:
    from distutils.core import setup


setup(**metadata)

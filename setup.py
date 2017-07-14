from __future__ import absolute_import
import re

versionstuff = dict(
    re.findall("(.+) = '(.+)'\n", open('rem3d/version.py').read()))

metadata = dict(name = 'rem3d',
                version=versionstuff['version'],
                description='a modeling and analysis toolkit for reference Earth datasets and tomographic models.',
                url='http://www.rem3d.org',
    		author = 'Pritwiraj Moulik',
                author_email='moulik.umd@gmail.com',
                license='GPL',
                long_description='REM3D is a Python library for reference Earth datasets and tomographic models.',
    			packages = ['rem3d'],
                package_data={'rem3d': ['data/input_*/*']},
    		keywords = ['earth-science','earth-observation','earthquake','earth','earthquake-data','geology','geophysics','geophysical-inversions','seismology','seismic-inversion','seismic-waves','seismic-tomography','mineral','geochemistry','geodesy','physics','modeling','modeling-tool','model','geodynamics'],
                classifiers=[
                'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                'Programming Language :: Python :: 2.7',
                'Programming Language :: Python :: 3.4'],
                )


'earth-science','earth-observation','earthquake','earth','earthquake-data','geology','geophysics','geophysical-inversions','seismology','seismic-inversion','seismic-waves','seismic-tomography','mineral','geochemistry','geodesy','physics','modeling','modeling-tool','model','geodynamics'
# Try to use setuptools in order to check dependencies.
# if the system does not have setuptools, fall back on
# distutils.
try:
    from setuptools import setup
    metadata['install_requires'] = ['numpy', 'matplotlib', 'scipy']
except ImportError:
    from distutils.core import setup


setup(**metadata)

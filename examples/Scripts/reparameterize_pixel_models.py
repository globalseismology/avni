#!/usr/bin/env python
"""This module contains the various subroutines used for plotting
Usage import """

import argparse #parsing arguments
import ntpath
import pdb
####################### IMPORT REM3D MODULES   #####################################
from rem3d.models.model3d import Model3D
#########################################################
def main():
    parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
    parser.add_argument('original', type=str,
        help='An orginial model file in REM3D netcdf format')
    parser.add_argument('target', type=str,
        help='A target model in REM3D netcdf format')
    parser.add_argument('-o', '--output', action='store_true',
        help='Save the figures in files')
    arg = parser.parse_args()

    # read the model to model3d instance
    print('.... reading original model from '+arg.original)
    original = Model3D(arg.original)
    print('.... reading target model from '+arg.target)
    target = Model3D(arg.target)

    # reparameterize the simulation to seismic
    reparam,tree = original.reparameterize(target,interpolated=True)

    # get the resolution matrix of the seismic from the same folder
    reparam.writerealization()
    return

if __name__== "__main__":
    main()

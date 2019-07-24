#!/usr/bin/env python
"""This module contains the various subroutines used for plotting
Usage import """

import argparse #parsing arguments
import ntpath
import pdb
####################### IMPORT REM3D MODULES   #####################################
from rem3d.tools import stage
from rem3d.data import update_file
from rem3d.models.model3d import Model3D
#########################################################
def main():
    parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
    parser.add_argument('simulation', type=str,
        help='A simulated geodynamic model file in REM3D netcdf format')
    parser.add_argument('seismic', type=str,
        help='A seismic model in REM3D netcdf format')
    parser.add_argument('-p', '--parameter', type=str, default='vs',
        help='Parameter of interest')
    parser.add_argument('-u', '--upper_bound', type=float, default=2.0,
        help='Upper bound for color scale saturation level (percent)')
    parser.add_argument('-l', '--lower_bound', type=float, default=-2.0,
        help='Lower bound for color scale saturation level (percent)')
    parser.add_argument('-i', '--elev_interval', type=float, default=100,
        help='Number of elevation points for transect plots')
    parser.add_argument('-e', '--elev_exxagerate', type=float, default=50,
        help='Elevation exxageration for transect plots')
    parser.add_argument('-t','--translate', help="translate physical parameters to seismic parameters", action="store_true")
    parser.add_argument('-o', '--output', action='store_true',
        help='Save the figures in files')
    arg = parser.parse_args()

    # read the model to model3d instance
    print('.... reading simulation from '+arg.simulation)
    simulation = Model3D(arg.simulation)
    print('.... reading seismic model from '+arg.seismic)
    seismic = Model3D(arg.seismic)
    pdb.set_trace()

    # reparameterize the simulation to seismic
    reparam = simulation.reparameterize(seismic)
    pdb.set_trace()

    # get the resolution matrix of the seismic from the same folder
    outsparse,modelarr = seismic.get_resolution()

    # m_filtered = R X m_reparam
    # m_out = outsparse * reparam

    # save and plot m_out reparam
    return

if __name__== "__main__":
    main()

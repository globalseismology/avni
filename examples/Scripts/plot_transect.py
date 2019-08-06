#!/usr/bin/env python
"""This module contains the various subroutines used for plotting
Usage import """

import argparse #parsing arguments
import ntpath
################################ IMPORT REM3D MODULES   #####################################
from rem3d.f2py import ddelazgc # geolib library from NSW
from rem3d.plots import plot1section
from rem3d.tools import stage
from rem3d.data import update_file
#########################################################
def main():
    parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
    parser.add_argument('-m', '--file', type=str, default='S362ANI+M.BOX25km_PIX1X1.rem3d.nc4',
        help='Model file')
    parser.add_argument('-p', '--parameter', type=str, default='vs',
        help='Parameter of interest')
    parser.add_argument('-u', '--upper_bound', type=float, default=2.0,
        help='Upper bound for color scale saturation level (percent)')
    parser.add_argument('-l', '--lower_bound', type=float, default=-2.0,
        help='Lower bound for color scale saturation level (percent)')
    parser.add_argument('-t', '--colorcontour', type=int, default=20,
        help='color contour levels')
    parser.add_argument('-i', '--elev_interval', type=float, default=100,
        help='Number of elevation points for transect plots')
    parser.add_argument('-e', '--elev_exxagerate', type=float, default=50,
        help='Elevation exxageration for transect plots')
    parser.add_argument('-o', '--output', action='store_true',
        help='Save the figures in files')
    parser.add_argument('-f', '--format', type=str, default='png',
        help='Outfile file format')
    arg = parser.parse_args()

    try:
        # stage the file for plotting
        ierror = stage(arg.file,overwrite=True)
        model3d = ntpath.basename(arg.file)
    except:
        # update the file from the server
        update_file(arg.file)
        model3d = arg.file

    ##### Example of a regional transects
    print("PLOTTING SECTION 1")
    # Kermadec
    lat1 = -25.;lng1 = 191.;lat2 = -22.;lng2 = 160.
    delta,azep, _  = ddelazgc(lat1,lng1,lat2,lng2)
    if arg.output:
        outfile = 'Kermadec.'+arg.format
    else:
        outfile = None
    topo,topo_tree,tomo,tomo_tree = plot1section(lat1,lng1,azep,delta,model=model3d,parameter=arg.parameter,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{S} / V_{S}$'+' (%)',outfile=outfile,vexaggerate=arg.elev_exxagerate,nelevinter=arg.elev_interval,colorcontour=arg.colorcontour)

    # N. Chile
    lat1 = -29.;lng1 = -50.;lat2 = -29.;lng2 = -80.
    delta,azep, _  = ddelazgc(lat1,lng1,lat2,lng2)
    if arg.output:
        outfile = 'NorthChile.'+arg.format
    else:
        outfile = None
    plot1section(lat1,lng1,azep,delta,topo=topo,topotree=topo_tree,modeltree=tomo_tree,model=tomo,parameter=arg.parameter,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{S} / V_{S}$'+' (%)',vexaggerate=arg.elev_exxagerate,nelevinter=arg.elev_interval,outfile=outfile,colorcontour=arg.colorcontour)

    #Japan
    lat1 = 34.;lng1 = 152.;lat2 = 40.;lng2 = 117.
    delta,azep,_ = ddelazgc(lat1,lng1,lat2,lng2)
    if arg.output:
        outfile = 'Japan.'+arg.format
    else:
        outfile = None
    plot1section(lat1,lng1,azep,delta,topo=topo,topotree=topo_tree,modeltree=tomo_tree,model=tomo,parameter=arg.parameter,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{S} / V_{S}$'+' (%)',vexaggerate=arg.elev_exxagerate,nelevinter=arg.elev_interval,outfile=outfile,colorcontour=arg.colorcontour)


    ###### Example of a 180 degree transect without topography
    print("PLOTTING SECTION 2")
    lat1 = 0.;lng1 = 0.;azimuth = -30.;gcdelta = 180.
    if arg.output:
        outfile = 'transect180.'+arg.format
    else:
        outfile = None
    plot1section(lat1,lng1,azimuth,gcdelta,topo=topo,topotree=topo_tree,modeltree=tomo_tree,model=tomo,parameter=arg.parameter,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{S} / V_{S}$'+' (%)',vexaggerate=0,figuresize=[8,4],width_ratios=[1,4],numevalx=360,numevalz=360,outfile=outfile,colorcontour=arg.colorcontour)

    ###### Example of a 360 degree transect without topography
    print("PLOTTING SECTION 3")
    lat1 = 0.;lng1 = 150.;azimuth = 90.;gcdelta = 360.
    if arg.output:
        outfile = 'transect360.'+arg.format
    else:
        outfile = None
    plot1section(lat1,lng1,azimuth,gcdelta,topo=topo,topotree=topo_tree,modeltree=tomo_tree,model=tomo,parameter=arg.parameter,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{S} / V_{S}$'+' (%)',vexaggerate=0,figuresize=[8,4],width_ratios=[1,4],numevalx=720,numevalz=720,outfile=outfile,colorcontour=arg.colorcontour)

# plot1section(lat1,lng1,azimuth,gcdelta,model=model3d,parameter=arg.parameter,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{S} / V_{S}$'+' (%)',vexaggerate=0,figuresize=[8,4],width_ratios=[1,4],numevalx=720,numevalz=720,colorcontour=arg.colorcontour)

    return

if __name__== "__main__":
    main()

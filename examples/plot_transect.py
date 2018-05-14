#!/usr/bin/env python
"""This module contains the various subroutines used for plotting
Usage import """

import sys,os
import argparse #parsing arguments

################################ IMPORT REM3D MODULES   #####################################
from rem3d.geolib import delazgc # geolib library from NSW
from rem3d.plots import plot1section
#########################################################
def main():
	parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
	parser.add_argument("file", type=str,
		help='Map file in .epix or 3D model file in .npz format')
	parser.add_argument('-u', '--upper_bound', type=float, default=1.0,
        help='Upper bound for color scale saturation level (percent)')
	parser.add_argument('-l', '--lower_bound', type=float, default=-1.0,
        help='Lower bound for color scale saturation level (percent)')
	parser.add_argument('-i', '--elev_interval', type=float, default=50,
        help='Number of elevation points for transect plots')
	parser.add_argument('-e', '--elev_exxagerate', type=float, default=0,
        help='Elevation exxageration for transect plots')
	parser.add_argument('-d', '--dbs_path', type=str, default='~/dbs',
		help='Database path containing files like hotspots.pkl')
	parser.add_argument('-b', '--base_path', type=str, default='.',
		help='Base path containing FILE')   
	parser.add_argument('-c', '--num_cores', type=int, default=1,
		help='Number of cores to use')
	arg = parser.parse_args()
	
	filename=arg.base_path+'/'+arg.file
# 	plot1globalmap(filename,arg.lower_bound,arg.upper_bound,arg.dbs_path,projection='robin',colorlabel='$\delta V_{P} / V_{P}$'+' (%)',lat_0=0,lon_0=150)
# 	plot1globalmap(filename,arg.lower_bound,arg.upper_bound,arg.dbs_path,projection='ortho',colorlabel='$\delta V_{P} / V_{P}$'+' (%)',lat_0=0,lon_0=190)
	
	
	##### Example of a regional transects
	#Japan
# 	lat1 = 34.;lng1 = 152.;lat2 = 40.;lng2 = 117.
# 	delta,azep,azst = delazgc(lat1,lng1,lat2,lng2)
# 	plot1section(lat1,lng1,azep,delta,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{P} / V_{P}$'+' (%)',vexaggerate=arg.elev_exxagerate,dbs_path=arg.dbs_path,width_ratios=[1,2],nelevinter=arg.elev_interval)
	
	# Kermadec
# 	lat1 = -25.;lng1 = 191.;lat2 = -22.;lng2 = 160.
# 	delta,azep,azst = delazgc(lat1,lng1,lat2,lng2)
# 	plot1section(lat1,lng1,azep,delta,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{P} / V_{P}$'+' (%)',vexaggerate=arg.elev_exxagerate,dbs_path=arg.dbs_path,width_ratios=[1,2],nelevinter=arg.elev_interval)

	# N. Chile	
	lat1 = -29.;lng1 = -50.;lat2 = -29.;lng2 = -80.
	delta,azep,azst = delazgc(lat1,lng1,lat2,lng2)
	plot1section(lat1,lng1,azep,delta,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{P} / V_{P}$'+' (%)',vexaggerate=50,dbs_path=arg.dbs_path,width_ratios=[1,2],nelevinter=arg.elev_interval)

	###### Example of a 180 degree transect without topography
	lat1 = 0.;lng1 = 0.;azimuth = -30.;gcdelta = 180.
	plot1section(lat1,lng1,azimuth,gcdelta,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{P} / V_{P}$'+' (%)',vexaggerate=arg.elev_exxagerate,dbs_path=arg.dbs_path,figuresize=[8,4],width_ratios=[1,4],nelevinter=arg.elev_interval)

	###### Example of a 360 degree transect without topography
	lat1 = 0.;lng1 = 0.;azimuth = -45.;gcdelta = 360.
	plot1section(lat1,lng1,azimuth,gcdelta,vmin=arg.lower_bound,vmax=arg.upper_bound,colorlabel='$\delta V_{P} / V_{P}$'+' (%)',vexaggerate=arg.elev_exxagerate,dbs_path=arg.dbs_path,figuresize=[8,4],width_ratios=[1,4],nelevinter=arg.elev_interval)
		
	return

if __name__== "__main__":
	main()

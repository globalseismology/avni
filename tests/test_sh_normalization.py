#!/usr/bin/env python
"""This module contains the various subroutines used for plotting
Usage import """

import argparse #parsing arguments
import ntpath
import numpy as np
import pyshtools
import matplotlib.pyplot as plt
import xarray as xr
import pdb

####################### IMPORT REM3D MODULES   #####################################

from rem3d.tools import wrswpsh,eval_ylm
####################### IMPORT REM3D MODULES   #####################################


def main():
    parser = argparse.ArgumentParser(description='compare spherical harmonic normalizations')
    parser.add_argument('-d','--harmonic', type=str, default='C2,0',
        help='Harmonic in Cl,m format')
    parser.add_argument('-s', '--spacing', type=float, default=5.0,
        help='Spacing in degrees')
    parser.add_argument('-o', '--output', action='store_true',
        help='Display figure')
    arg = parser.parse_args()

    spacing = arg.spacing
    harmonic = arg.harmonic

    data_vars={}
    L, M = harmonic[1:].split(',')
    CS = harmonic[0]; L = int(L); M = int(M)
    if CS.lower() == 's' and M==0: raise ValueError('no sine value for order 0')
    lmaxhor = L
    # find index
    loop = np.power(L,2) if L > 0 else 0
    index = None
    if L == 0:
        index =loop
    else:
        for order in np.arange(L+1):
            cslist = ['C'] if order == 0 else ['C','S']
            for cs in cslist:
                if order == M and CS.lower() == cs.lower():
                    index = loop
                    break
                loop = loop + 1

    # In pyshtools:
    # A 2D equally sampled (default) or equally spaced map in degrees of the spherical
    # harmonic coefficients cilm that conforms to the sampling theorem of Driscoll and
    # Healy (1994). The first latitudinal band corresponds to 90 N, the latitudinal band
    # for 90 S is not included, and the latitudinal sampling interval is 180/n degrees,
    # where n is 2*lmax+2. The first longitudinal band is 0 E, the longitudinal band for
    # 360 E is not included, and the longitudinal sampling interval is 360/n for an
    # equally sampled and 180/n for an equally spaced grid, respectively.
    latitude = np.arange(90.,-90.,-spacing)
    longitude = np.arange(0.,360.,spacing)

    # ylm values
    ylmcof = eval_ylm(latitude,longitude,lmaxhor,grid=True)
    data_vars['ylm'] = (('latitude', 'longitude'),(ylmcof[:,index].reshape(len(latitude),len(longitude))).toarray())

    # shold values with coefficients scaled to account for normalization
    sholdcof = eval_ylm(latitude,longitude,lmaxhor,grid=True,norm='shold')
    if M > 0: sholdcof = sholdcof/np.sqrt(2.)
    data_vars['shold'] = (('latitude', 'longitude'),(sholdcof[:,index].reshape(len(latitude),len(longitude))).toarray())

    # pyshtools with coefficients scaled to account for normalization
    cilm = np.zeros((2,lmaxhor+1,lmaxhor+1))
    if CS.lower()=='c' and M == 0: cilm[0,L,M]=1./np.sqrt(4*np.pi)
    if M > 0:
        if CS.lower()=='c' :cilm[0,L,M]=np.power(-1.,M)/np.sqrt(8*np.pi)
        if CS.lower()=='s' : cilm[1,L,M]=np.power(-1.,M)/np.sqrt(8*np.pi)
    topo = pyshtools.expand.MakeGridDH(cilm, sampling=2,norm=1,lmax=(180./spacing-2)/2)
    data_vars['shtools'] = (('latitude', 'longitude'),topo)
    xrdata = xr.Dataset(data_vars=data_vars,coords = {'latitude':latitude,'longitude':longitude})
    # checks
    expected = xrdata.dims['latitude'] * xrdata.dims['longitude']
    diff1 = xrdata['shtools']-xrdata['ylm']
    diff2 = xrdata['shold']-xrdata['ylm']
    assert np.sum(np.isclose(diff1,0,atol=1e-5)) == expected,' values of shtools and ylm should be similar'
    assert np.sum(np.isclose(diff2,0,atol=1e-5)) == expected,' values of shold and ylm should be similar'
    if arg.output:
        fig, axes = plt.subplots(nrows=4,sharex='col')
        xrdata['ylm'].plot(ax=axes[0])
        xrdata['shold'].plot(ax=axes[1])
        xrdata['shtools'].plot(ax=axes[2])
        (diff1).plot(ax=axes[3])
        axes[3].set_title('Difference (shtools-ylm)')
        fig.suptitle('Spherical harmonic: '+harmonic, fontsize=16)
        plt.show()

if __name__== "__main__":
    main()

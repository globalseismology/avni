#!/usr/bin/env python

#####################  IMPORT STANDARD MODULES   #########################

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import float,int,list,tuple

import os
import numpy as np
from math import atan
from configobj import ConfigObj
import typing as tp

####################### IMPORT AVNI LIBRARIES  ###########################

from .. import constants
from ..tools import convert2units,get_configdir,get_fullpath

##########################################################################

def getplanetconstants(planet: tp.Union[None,str] = None, configfile: tp.Union[None,str] = None, option = None):
    """Load the astronomic-geodetic constraints for a planet from a
    configuration file.

    Parameters
    ----------
    planet : tp.Union[None,str], optional
        _description_, by default None so :py:func:`constants.planetpreferred`
    configfile : tp.Union[None,str], optional
        all the planet configurations are in this file., by default None so
        read from so :py:func:`get_configdir/constants.planetpreferred`
    option
        GRS option for constants, by default None so use the default one

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # defaults
    if planet is None: planet = constants.planetpreferred
    if configfile is None: configfile = get_fullpath(get_configdir()+'/'+constants.planetconstants)

    if not os.path.isfile(configfile):
        raise IOError('No configuration file found: '+configfile)
    else:
        parser = ConfigObj(configfile)

    if option == None:
        try:
            option = parser[planet]['default']
        except:
            raise KeyError('Default GRS needs to be specified for '+planet+' in file '+configfile+'. Or an option passed as an argument to getplanetconstants')

    try:
        parser_select = parser[planet][option]
    except:
        raise IOError('No planet '+planet+' with GRS '+option+' found in file '+configfile)
    constants.a_e = convert2units(parser_select['a_e']) # Equatorial radius
    constants.GM = convert2units(parser_select['GM']) # Geocentric gravitational constant m^3s^-2
    constants.G = convert2units(parser_select['G']) # Gravitational constant m^3kg^-1s^-2
    constants.R = convert2units(parser_select['R']) # Radius of the Earth in m
    try:
        constants.f = convert2units(parser_select['f']) #flattening
        constants.rf = 1./constants.f #inverse flattening
    except KeyError:
        try:
            constants.rf = convert2units(parser_select['1/f']) #inverse flattening
            constants.f = 1./constants.rf #flattening
        except:
            raise KeyError('need either flattening (f) or inverse flattening (1/f) for '+planet+' in '+configfile)
    constants.omega = convert2units(parser_select['omega']) #Angular velocity in rad/s
    constants.M_true = convert2units(parser_select['M_true']) # Solid Earth mass in kg
    constants.I_true = convert2units(parser_select['I_true'])# Moment of inertia in m^2 kg
    constants.rhobar = convert2units(parser_select['rhobar']) # Average density in kg/m^3
    #length of 1 degree in km
    try:
        constants.deg2km = convert2units(parser_select['deg2km'])
    except KeyError:
        constants.deg2km = 2*np.pi*constants.a_e/360/1000
    constants.deg2m = constants.deg2km * 1000. #length of 1 degree in m
    #DeltaJ2 Permanent tide correction
    try:
        constants.DeltaJ2 = convert2units(parser_select['DeltaJ2'])
    except KeyError:
        constants.DeltaJ2 = None
    try:
        constants.k2 = convert2units(parser_select['k2'])
    except KeyError:
        constants.k2 = None
    try:
        constants.barC2hydro = convert2units(parser_select['barC2hydro'])
    except KeyError:
        constants.barC2hydro = None
    try:
        constants.barC4hydro = convert2units(parser_select['barC4hydro'])
    except KeyError:
        constants.barC4hydro = None
    # correction for geographic-geocentric conversion: 0.993277 for 1/f=297
    try:
        print('... Re - Initialized avni module with constants for '+planet+' from '+parser_select['cite']+' from geocentric correction '+str(constants.geoco))
        constants.geoco = (1.0 - constants.f)**2.
    except AttributeError:
        constants.geoco = (1.0 - constants.f)**2.
    constants.planet = planet
    constants.grs = option


def evaluate_grs(GM: tp.Union[None,float] = None,
                 f: tp.Union[None,float] = None,
                 a_e: tp.Union[None,float] = None,
                 omega: tp.Union[None,float] = None,
                 R: tp.Union[None,float] = None,
                 nzo: int = 10,
                 store: bool = False):
    """Calculate geopotential constants in a reference earth model.

    All the following page numbers and equation numbers refer to the
    book Physical Geodesy by Hofmann-wellenhof and Moritz :cite:p:`hofmann2006physical`

    Parameters
    ----------
    GM : tp.Union[None,float], optional
        Gravitational constant times mass reference, by default None
    f : tp.Union[None,float], optional
        Flattening, by default None
    a_e : tp.Union[None,float], optional
        Semi-major axis, by default None
    omega : tp.Union[None,float], optional
        Angular velocity, by default None
    R : tp.Union[None,float], optional
        _description_, by default None
    nzo : int, optional
        Number of zonal harmonics (2,4,... 2*nzo), by default 10
    store : bool, optional
        Store in constants or return as output if False, by default False

    Returns
    -------
    barC2n
        Normalized even zonal harmonics of
        the corresponding Somigliana-Pizzetti normal field.
        barC2n(:,1): normalized zonal harmonics
        barC2n(:,2): degree of the zonal harmonic [2 4 ... 2*nzo]
    geqt
        Normal gravity at the equator
    gpol
        Normal gravity at the pole
    U0
        Normal potential at the ellipsoid
    m
        omega^2*a^2*b/(GM)
    ecc
        First eccentricity
    eccp
        Second eccentricity
    a_p
        Semi-minor axis
    E
        Linear eccentricity
    c
        Polar radius of curvature

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # defaults are from constants
    if GM == None: GM = constants.GM
    if f == None: f = constants.f
    if a_e == None: a_e = constants.a_e
    if R == None: R = constants.R
    if omega == None: omega = constants.omega

    # First eccentricity;
    ecc=np.sqrt(2*f-f**2)
    if store: constants.ecc = ecc

    # First eccentricity squared
    ecc2=ecc**2
    if store: constants.ecc2 = ecc2

    # Second eccentricity
    # p. 71, Eqn.(2-138);
    eccp=ecc/(1-f)
    if store: constants.eccp = eccp
    # Second eccentricity squared
    eccp2=eccp**2
    if store: constants.eccp2 = eccp2

    # Semi-minor axis
    a_p=a_e*(1-f)
    if store: constants.a_p = a_p

    # m
    # p. 70, Eqn.(2-137);
    m=(omega**2)*(a_e**2)*a_p/GM;
    if store: constants.m = m

    # q_0 and q_0p
    # p. 67, Eqn.(2-113)
    # In the book, q_0 is q
    q_0=1/2*((1+3/eccp2)*atan(eccp)-3/eccp)
    if store: constants.q_0 = q_0
    # and q_0p is q_0
    q_0p=3*(1+1/eccp2)*(1-1/eccp*atan(eccp))-1
    if store: constants.q_0p = q_0p

    # J_2
    # p. 75-76, Eqn.(2-165), Eqn.(2-166) and Eqn.(2-172)
    j_2=ecc2/3*(1-2/15*m*eccp/q_0)
    if store: constants.j_2 = j_2

    # j_2n
    # p. 76, Eqn.(2-170) and Eqn.(2-172)
    j_2n = {}
    for N in np.arange(1,nzo+1):
        j_2n[2*N]=((-1)**(N+1))*3*(ecc2**N)/(2*N+1)/(2*N+3)*(1-N+5*N*j_2/ecc2)
    if store: constants.j_2n = j_2n

    # Normalized C2n0 terms.
    # p. 60, Eqn.(2-80)
    barC2n={}
    for N in np.arange(1,nzo+1): barC2n[2*N] = -j_2n[2*N]/np.sqrt(2*2*N+1)
    if store: constants.barC2n = barC2n

    # Normal gravity at the equator.
    # p. 71, Eqn.(2-141);
    geqt=GM/a_e/a_p*(1-m-m/6*eccp*q_0p/q_0)
    if store: constants.geqt = geqt

    # Normal gravity at the pole.
    # p. 71, Eqn.(2-142)
    gpol=GM/(a_e**2)*(1+m/3*eccp*q_0p/q_0)
    if store: constants.gpol = gpol

    # Mean Normal gravity Somiglinana's closed formula at phi=45.
    k = (a_p*gpol/a_e/geqt)-1
    g0= geqt*(1+k/2)/np.sqrt(1-ecc2/2)
    if store: constants.g0 = g0

    # Mean gravitational acceleration on the sphere
    g=GM/(R**2)
    if store: constants.g = g

    # Linear eccentricity
    # p. 66, Eqn.(2-101)
    E=np.sqrt(a_e**2-a_p**2);
    if store: constants.E = E

    # Normal potential at the ellipsoid
    # p. 68, Eqn.(2-123)
    U0=GM/E*atan(eccp)+1/3*(omega**2)*(a_e**2)
    if store: constants.U0 = U0

    # Polar radius of curvature
    # p. 73, eqn.(2-150)
    c=(a_e**2)/a_p
    if store: constants.c = c

    if not store: return barC2n,j_2n,geqt,gpol,g0,g,U0,m,ecc,eccp,a_p,E,c
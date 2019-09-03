#!/usr/bin/env python
"""
This module contains the various subroutines used for plotting
Usage import
"""
#####################  IMPORT STANDARD MODULES   ######################################

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

import os
import numpy as np #for numerical analysis
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import colorsys
import traceback
import pdb

####################       IMPORT OWN MODULES     ######################################
from .. import tools
from .. import data
from .. import constants
########################      GENERIC   ################################################

def updatefont(fontsize=15,fontname='sans-serif',ax=None):
    """
    Updates the font type and sizes globally or for a particular axis handle

    Parameters
    ----------

    ax :  figure axis handle

    fontsize,fontname : font parameters

    Return:
    ----------

    ax : updated axis handle if ax is not None

    """
    if ax is None:
        plt.rcParams["font.family"] = fontname
        plt.rcParams["font.size"] = fontsize
    else:
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
            item.set_fontname(fontname)
        for item in ([ax.xaxis.label, ax.yaxis.label]):
            item.set_fontsize(fontsize+2)
            item.set_fontname(fontname)
        ax.title.set_fontsize(fontsize+3)
        ax.title.set_fontname(fontname)
    return ax if ax is not None else None

def initializecolor(name,**kwargs):
    """
    Initialize a color palette instance from standard Python, constants.py or downloadable from server

    name : name of color palette.Can have _r appended to standard ones
                    for reversed color scales

    kwargs : optional arguments for Basemap
    """
    success1 = True; success2 = True; success3 = True # assume colorscale is read somehow
    try:
        cpalette = plt.get_cmap(name)
    except ValueError:
        success1 = False
        var1 = traceback.format_exc()
        try: # try getting color palette from standard ones in constants.py
            cpalette = standardcolorpalette(name)
        except KeyError:
            success2 = False
            var2 = traceback.format_exc()
            try:
                if kwargs:
                    cpalette = customcolorpalette(name,**kwargs)
                else:
                    cpalette = customcolorpalette(name)
            except:
                success3 = False
                print('#########   Tried reading as standard Python palette   ##########')
                print(var1)
                print('#########   Tried reading as  ones from constant.py   ##########')
                print(var2)
                print('############    Tried downloading from server   ############')
                print(traceback.format_exc())
    if not success1 and not success2 and not success3:
        raise IOError('unable to read color palette '+name+' from standard Python, constants.py or downloadable from server.')
    return cpalette

def standardcolorpalette(name='rem3d'):
    """
    Get a custom REM3D color palette from constants.py

    Parameters
    ----------

    name : color palette name that will be used elsewhere
           if name ends in '_r', use the reversed color scale.
    reverse: if the colors need to be reversed from those provided.

    """
    if name.endswith('_r'):
        RGBoption = name.split('_r')[0]
        RGBlist=constants.colorscale[RGBoption]['RGB'][::-1]
    else:
        RGBoption = name
        RGBlist=constants.colorscale[RGBoption]['RGB']
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(name, RGBlist,N=len(RGBlist))
    cmx.register_cmap(name=custom_cmap.name, cmap=custom_cmap)
    return custom_cmap

def get_colors(val,xmin=-1.,xmax=1.,palette='coolwarm',colorcontour=20):
    """gets the value of color for a given palette"""
    cm = cmx.get_cmap(palette)
    #cNorm  = mcolors.Normalize(vmin=xmin, vmax=xmax)
    bounds = np.linspace(xmin,xmax,colorcontour+1)
    cNorm = mcolors.BoundaryNorm(bounds,cm.N)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=palette)
    colorVal = scalarMap.to_rgba(val)
    return colorVal

def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))

    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]

    return cmap.from_list(cmap.name + "_gray", colors, cmap.N)

def make_colormap(seq,name='CustomMap'):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap(name, cdict)

def getcolorlist(cptfile,type='rem'):
    """Get a tuple for colorlist from a cptfile"""

    if not os.path.isfile(cptfile): raise IOError("File ("+cptfile+") does not exist.")
    colorlist=[]

    # if it is the format in the REM project
    if type=='rem':
        cptarr=np.genfromtxt(cptfile, dtype=None,comments="#")
        for irow in np.arange(len(cptarr)):
            tups=cptarr[irow][1]/255.,cptarr[irow][2]/255.,cptarr[irow][3]/255.
            val=(cptarr[irow][4]-cptarr[0][4])/(cptarr[len(cptarr)-1][0]-cptarr[0][4])
            if irow==1:
                colorlist.append(tups)
            elif irow > 1 and irow < len(cptarr)-1:
                colorlist.append(tups)
                colorlist.append(val)
                colorlist.append(tups)

    # if it is the standard format in the GMT project
    elif type=='standard':
        colorlist = readstandardcpt(cptfile)

    return colorlist

def readstandardcpt(cptfile):
    """Read a GMT color map from an OPEN cpt file
    Parameters
    ----------
    cptf : open file or url handle
        path to .cpt file
    name : str, optional
        name for color map
        if not provided, the file name will be used
    """
    if not os.path.isfile(cptfile): raise IOError("File ("+cptfile+") does not exist.")

    # process file
    x = []; r = []; g = []; b = []
    lastls = None

    fo = open(cptfile, "r")
    cptlines = fo.readlines()
    colorlist = []

    for l in cptlines:
        ls = l.split()

        # skip empty lines
        if not ls: continue

        # parse header info
        if ls[0] in ["#", b"#"]:
            if ls[-1] in ["HSV", b"HSV"]:
                colorModel = "HSV"
            else:
                colorModel = "RGB"
            continue

        # skip BFN info
        if ls[0] in ["B", b"B", "F", b"F", "N", b"N"]: continue

        # parse color vectors
        x.append(float(ls[0]))
        r.append(float(ls[1]))
        g.append(float(ls[2]))
        b.append(float(ls[3]))

        # save last row
        lastls = ls

    x.append(float(lastls[4]))
    r.append(float(lastls[5]))
    g.append(float(lastls[6]))
    b.append(float(lastls[7]))

    if colorModel == "HSV":
        for i in range(len(r)):
            # convert HSV to RGB
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360., g[i], b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    elif colorModel == "RGB":
        r = [val/255. for val in r]
        g = [val/255. for val in g]
        b = [val/255. for val in b]

    x = np.array(x)
    xNorm = (x - x[0])/(x[-1] - x[0])
    for i in range(len(x)):
        tups=r[i],g[i],b[i]
        colorlist.append(tups)
        colorlist.append(xNorm[i])
        colorlist.append(tups)

    # Close opend file
    fo.close()

    # return colormap
    return colorlist

def downloadcpt(url):
    file = url.split('/')[-1]
    baseurl = url.replace('/'+file,'')
    data.update_file(file,folder=tools.get_cptdir(),baseurl = baseurl)

def customcolorpalette(name='bk',cptfolder=tools.get_cptdir(),colorlist=None,colormax=2.,middlelimit=0.5,ifgraytest=0):
    """Used to return preset color palettes from cptfolder. ifgraytest test how the figure looks in gray scale. (-colormax,colormax) are the limits of the colorbar. zerolimit is the limit to which the middle color (e.g. grey) will extend on either side of colorttmax mid. """
    c = mcolors.ColorConverter().to_rgb
    if name=='r_lgrey_b':
        colorlist=[c('blue'), c('lightgray'), (2.*colormax-2.*middlelimit)/(4.*colormax), c('lightgray'),c('lightgray'), (2.*colormax+2.*middlelimit)/(4.*colormax), c('lightgray'),c('red'), 1., c('red')]
    else:
        if name=='bk':
            file = 'bk1_0.cpt_'
            type = 'rem'
        elif name=='hit1':
            file = 'hit1.cpt_'
            type = 'rem'
        elif name=='yuguinv':
            file = 'yu1_2inv.new.cpt_'
            type = 'rem'
        else:
            file = name+'.cpt'
            type = 'standard'
        if type =='rem': data.update_file(file,folder=cptfolder)
        cptfolder=tools.get_fullpath(cptfolder)
        if os.path.isfile(cptfolder+'/'+file):
            colorlist=getcolorlist(cptfolder+'/'+file,type=type)
        else:
            raise ValueError("Could not find file "+cptfolder+'/'+file)

    if colorlist is None: raise ValueError("No colorlist found")
    custom_cmap = make_colormap(colorlist,name)
    cmx.register_cmap(name=custom_cmap.name, cmap=custom_cmap)
    palette=custom_cmap.name

    if ifgraytest==1:
        palette=grayify_cmap(palette)

    return custom_cmap

from __future__ import absolute_import
import scipy.constants
import ConfigParser
import pdb    #for the debugger pdb.set_trace()
import pkgutil
import os
import sys
import codecs,json #printing output
import numpy as np
import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import pandas as pd
from math import ceil
####################### IMPORT REM3D LIBRARIES  #######################################
from . import constants
#######################################################################################

def get_fullpath(path):
    """Provides the full path by replacing . and ~ in path."""
    # Get the current directory    
    if path[0]=='.': path = os.path.dirname(os.path.abspath(__file__))+path[1:]   
    # If the path starts with tilde, replace with home directory
    if path[0]=='~': path=os.path.expanduser("~")+path[1:]
    return path
    
def listfolders(path):
    """
    Return a list of directories in a path
    """
    dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return dirs

def get_installdir(checkwrite=True,checkenv=True):
    """
    Get the installation directory for the rem3d module. checkwrite checks for write access to the files.
    checkenv checks if the directory is specified as an environment variable.
    """
    installdir='not_set'
    ierror=0
    if checkenv:
        if os.environ.get('REM3Ddir') is not None:
            installdir=os.environ.get('REM3Ddir')
            print "Warning: Reading REM3Ddir from environment variables. "

    if installdir == 'not_set':
        loader=pkgutil.find_loader('rem3d')
        if loader is None:
            print "Warning: installation directory not found for rem3d. Using current directory "
            installdir = os.getcwd()
        else:
            installdir = os.path.dirname(loader.get_filename('rem3d'))
            if installdir is None:
                print "Error: Specify REM3Ddir environment variable (export REM3Ddir=/path/to/rem3d/) or install rem3d module. "
                ierror=1
            if checkwrite:
                if not os.access(installdir, os.W_OK):
                    print "Error: Cannot I/O to rem3d directory "+installdir+ " due to permissions. \
Specify rem3d_dir environment variable or chdir to a different directory with I/O access. "
                    ierror=2
    return installdir, ierror

def get_filedir(checkwrite=True,makedir=True):
    """
    Get the local files directory. Make a new directory if doesn't exist (makedir==True)
    """
    installdir,ierror = get_installdir(checkwrite=checkwrite)
    filedir = installdir+'/'+constants.localfilefolder
    if checkwrite and makedir: 
        if not os.path.exists(filedir):
            os.makedirs(filedir)        
    return filedir
    
def writejson(nparray,filename,encoding='utf-8'):
    """Writes a json file from a numpy array"""
    
    listarray = nparray.tolist() # nested lists with same data, indices
    json.dump(listarray, codecs.open(filename, 'w', encoding=encoding), separators=(',', ':'), sort_keys=True, indent=4) ### this saves the array in .json format
    return

def readjson(filename,encoding='utf-8'):
    """Reading from a filename to a numpy array"""
        
    obj_text = codecs.open(filename, 'r', encoding=encoding).read()
    listarray = json.loads(obj_text)
    nparray = np.array(listarray)
    return nparray    

def uniquenumpyrow(a):
    """Gets the unique rows from a numpy array and the indices. e.g. to get unique lat-lon values"""
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    unique_a = a[idx]
    return unique_a,idx
<<<<<<< HEAD
=======
    
def interp_weights(xyz, uvw):
    """
    Gets the interpolation weights for a grid
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
    """
    tri = qhull.Delaunay(xyz)
    simplex = tri.find_simplex(uvw)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uvw - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts):
    """
    Makes the interpolation weights for a previous grid weights
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
    """
    return np.einsum('nj,nj->n', np.take(values, vtx), wts)
>>>>>>> 4c0183250a7f880765e859c035eb0b05134b4378

def getcommonSWcatalogs(SWdata1,SWdata2,shortref1='SW1',shortref2='SW2',decimals=2,write_common=False):
    """Get common path data for scatter plot of two SW catlog data. Lat/lon are rounded to decimal points."""
    patharr1=np.vstack((SWdata1['stlat'],SWdata1['stlon'],SWdata1['eplat'],SWdata1['eplon'])).T
    patharr2=np.vstack((SWdata2['stlat'],SWdata2['stlon'],SWdata2['eplat'],SWdata2['eplon'])).T

    # Round off to one decimal point to aid comparison
    patharr1=np.around(patharr1,decimals=decimals)
    patharr2=np.around(patharr2,decimals=decimals)

    # Get the unique combinations
    uni1,_ =uniquenumpyrow(patharr1)
    uni2,_ =uniquenumpyrow(patharr2)

    # Get the common unique paths between two catalogs
    nrows, ncols = uni1.shape
    dtype={'names':['f{}'.format(i) for i in range(ncols)],'formats':ncols * [uni1.dtype]}
    C = np.intersect1d(uni1.view(dtype), uni2.view(dtype))
    # This last bit is optional if you're okay with "C" being a structured array.
    C = C.view(uni1.dtype).reshape(-1, ncols)

    # now loop over common paths and take out delphase and distkm
    obs1=[None]*C.shape[0];obs2=[None]*C.shape[0]
    dist1=[None]*C.shape[0];dist2=[None]*C.shape[0]
    reportphase1=[None]*C.shape[0];reportphase2=[None]*C.shape[0]
    reporterr1=[None]*C.shape[0];reporterr2=[None]*C.shape[0]
    stat=[None]*C.shape[0];stlon=[None]*C.shape[0];stlat=[None]*C.shape[0]
    eplon=[None]*C.shape[0];eplat=[None]*C.shape[0]

    numcommon=C.shape[0]
    print '...... '+str(numcommon)+' measurements in common'
    if C.shape[0]>0:
        #bar = progressbar.ProgressBar(maxval=C.shape[0], \
        #widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        #bar.start()
        for ii in np.arange(C.shape[0]):
            test11=patharr1[:,0]==C[ii][0]
            test12=patharr1[:,1]==C[ii][1]
            test13=patharr1[:,2]==C[ii][2]
            test14=patharr1[:,3]==C[ii][3]
            indfind=np.where(np.logical_and(np.logical_and(np.logical_and(test11,test12),test13),test14))[0]
            obs1[ii]=SWdata1['delobsphase'][indfind].tolist()[0]
            dist1[ii]=SWdata1['distkm'][indfind].tolist()[0]
            reportphase1[ii]=SWdata1['refphase'][indfind].tolist()[0]
            reporterr1[ii]=SWdata1['delerrphase'][indfind].tolist()[0]
            stat[ii]=SWdata1['stat'][indfind].tolist()[0]
            stlat[ii]=SWdata1['stlat'][indfind].tolist()[0]
            stlon[ii]=SWdata1['stlon'][indfind].tolist()[0]
            eplat[ii]=SWdata1['eplat'][indfind].tolist()[0]
            eplon[ii]=SWdata1['eplon'][indfind].tolist()[0]
            #
            test21=patharr2[:,0]==C[ii][0]
            test22=patharr2[:,1]==C[ii][1]
            test23=patharr2[:,2]==C[ii][2]
            test24=patharr2[:,3]==C[ii][3]
            indfind=np.where(np.logical_and(np.logical_and(np.logical_and(test21,test22),test23),test24))[0]
            obs2[ii]=SWdata2['delobsphase'][indfind].tolist()[0]
            dist2[ii]=SWdata2['distkm'][indfind].tolist()[0]
            reportphase2[ii]=SWdata2['refphase'][indfind].tolist()[0]
            reporterr2[ii]=SWdata2['delerrphase'][indfind].tolist()[0]
            #bar.update(ii+1)
        #bar.finish()

        if write_common:
            filename=shortref1+'.'+shortref2+'.'+str(SWdata1['overtone'][0])+'.'+str(int(SWdata1['peri'][0]))+'.'+SWdata1['typeiorb'][0]+'.save'
            f = open(filename, 'wb')
            for obj in [stat,stlat,stlon,eplat,eplon,obs1, obs2, dist1, dist2, reportphase1, reportphase2, reporterr1, reporterr2]:
               cPickle.dump(obj, f, protocol=cPickle.HIGHEST_PROTOCOL)
            f.close()

    return numcommon,obs1,obs2,dist1,dist2

def getcommonSWcatalogsPandas(SWData1,SWData2,DegSig=2):
    # pulls common eq-station paths from two surface wave data sets

    # round the lat/lon columns for matching
    SWData1=roundSWData(SWData1,DegSig)
    SWData2=roundSWData(SWData2,DegSig)

    # merge the two dataframes, drop duplicates to get unique station-eq pairs
    col_list=['stlat','stlon','eplat','eplon']
    xtra_cols=['cmtdep','stat','distkm','delobsphase']

    # pd will append _x and _y to columns after merge if they exist in both df's
    xtra_x = [x + '_1' for x in xtra_cols]
    xtra_y = [x + '_2' for x in xtra_cols]
    xtra_cols = xtra_x+xtra_y
    xtra_cols.sort()

    # join them!
    SW_StatEpPairs=pd.merge(SWData1,SWData2,how='inner',on=col_list,suffixes=('_1', '_2'))
    SW_StatEpPairs=SW_StatEpPairs[col_list+xtra_cols] # drop other columns
    print('Records matching between groups: '+str(len(SW_StatEpPairs)))

    # get unique paths, stations and epicenters
    unique_stats=SW_StatEpPairs.copy(deep=True)
    unique_stats.drop_duplicates(subset=['stlat','stlon'],inplace=True)
    unique_stats['type']='station'
    print('Unique stations matching between groups: '+str(len(unique_stats)))
    unique_eps=SW_StatEpPairs.copy(deep=True)
    unique_eps.drop_duplicates(inplace=True,subset=['eplat','eplon'])
    unique_eps['type']='epicenter'
    print('Unique epicenters matching between groups: '+str(len(unique_eps)))

    return SW_StatEpPairs,unique_stats,unique_eps

def compareSWStudiesPandas(SW_StatEpPairs,deg_bin_width=2):

    # calculate comparisons between paths
    d2deg=1.0 / 111.1949
    SW_StatEpPairs['distdeg']=SW_StatEpPairs.distkm_1*d2deg
    SW_StatEpPairs['absdiff']=abs(SW_StatEpPairs['delobsphase_2']-SW_StatEpPairs['delobsphase_1'])

    # define bins
    nbins=ceil((SW_StatEpPairs.distdeg.max()-SW_StatEpPairs.distdeg.min())/deg_bin_width)
    bins=np.linspace(SW_StatEpPairs.distdeg.min(),SW_StatEpPairs.distdeg.max(),nbins)

    # calculate stats on binned values
    aggdict={'absdiff':['mean','median','std','count']}
    rndict={'mean':'mean_','median':'median_','std':'std_','count':'count_'}
    grpd = (SW_StatEpPairs.groupby
            (pd.cut(SW_StatEpPairs.distdeg,bins),as_index=False)
            .agg(aggdict)
            .rename(rndict))
    grpd['bin_center']=(bins[0:-1]+bins[1:])/2

    return SW_StatEpPairs,grpd
  
def roundSWData(SWdf,deg=2):
    # deg: degree decimal to round to (e.g., 1 = 0.1 degree bin), could be input
    round_dict={'stlat':deg,'stlon':deg,'eplat':deg,'eplon':deg}
    return SWdf.round(round_dict)

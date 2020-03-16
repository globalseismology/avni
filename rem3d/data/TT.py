#!/usr/bin/env python

#######################################################################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

#####################  IMPORT STANDARD MODULES   ######################################

import os
import sys
import h5py
import glob
import numpy as np
import pandas as pd
from itertools import islice

if sys.version_info[0] >= 3: unicode = str

####################       I/O ROUTINES     ######################################

def writetablehdf5(cagc_rays_outdir,tt_table_name,source_depth,component,verbose=True):
    '''
    writes an hdf5 file based on the output of CAGCRAYS

    params:
    cagc_rays <str>: path to output directory of CAGCRAYS
    tt_table_name <str>: name of hdf5 file to create or add to
    source_depth <float>: source depth in km
    component <str>: 'PSV', or 'SH'
    verbose <bool>: print extra output, helpful for debugging
    '''

    bad_paths = open('bad_paths.txt','w')
    source_depth = str(source_depth)

    #open hdf5 file
    if os.path.exists(tt_table_name):
        ds = h5py.File(tt_table_name,'r+')
    else:
        ds = h5py.File(tt_table_name,'w')

        #create h5py groups for travel times and paths
        ds.create_group('travel_times')
        ds['travel_times'].create_group('1D')
        ds.create_group('paths')

    try:
        ds['travel_times']['1D'].create_group(component)
    except:
        print('group already exists... appending new data to it')


    #get list of directories containing travel time and path into
    phase_list = glob.glob(cagc_rays_outdir+'/*')
    print(phase_list)


    for phase_ in phase_list:
        phase = phase_.split('/')[-1]
        print(phase)

        #get travel times
        try:
            tts = pd.read_csv(phase_+'/ttandtstar.txt',
                              delim_whitespace=True,
                              error_bad_lines=False,
                              skiprows=1,
                              names=['delta','p','T','dDdp','tstar','branch'])
            tts = tts.iloc[::-1]
            tts = tts.dropna() #get rid of NaN values

        except(OSError):
            raise ValueError('The directory ',phase,'does not contain ttandtstar.txt')

        #create a new group for each phase
        try:
            ds['travel_times']['1D'][component].create_group(phase)
        except:
            print('group already exists... appending new data to it')

        try:
            ds['travel_times']['1D'][component][phase].create_group(source_depth)
        except:
            print('group already exists... appending new data to it')

        #create datasets
        ds['travel_times']['1D'][component][phase][source_depth].create_dataset('distance_in_degree',data=tts['delta'].astype(float))
        ds['travel_times']['1D'][component][phase][source_depth].create_dataset('ray_param_sec_degree',data=tts['p'].astype(float))
        ds['travel_times']['1D'][component][phase][source_depth].create_dataset('time',data=tts['T'].astype(float))
        ds['travel_times']['1D'][component][phase][source_depth].create_dataset('dDdp',data=tts['dDdp'].astype(float))
        ds['travel_times']['1D'][component][phase][source_depth].create_dataset('tstar',data=tts['tstar'].astype(float))
        ds['travel_times']['1D'][component][phase][source_depth].create_dataset('branch',data=tts['branch'].values.astype('|S2'))

        #get raypaths
        path_list = glob.glob(cagc_rays_outdir+'/'+phase+'/path_SYNTH*')
        path_list.sort()

        try:
            ds['paths'].create_group(phase)
        except:
            print('group already exists... appending new data to it')

        for path_ in path_list:

            #read header
            with open(path_) as file_:
                header = list(islice(file_, 8))
            path_info = header[6]
            gcarc = path_info.strip().split()[1]

            #read data
            print("#########################################################")
            print(path_)
            path = np.loadtxt(path_)

            #get branch info
            with open(path_) as file_:
                branches = []
                for line_ in file_:
                    if line_.startswith("#BRANCH"):
                        branches.append(line_.strip().split()[1])

            #skip if there are two of the same branch at a single distance
            if len(branches) != len(set(branches)):
                bad_paths.write('source depth {}, phase {}, distance {}\n'.format(source_depth,phase,gcarc))
                continue

            if verbose:
                print('THE BRANCHES ARE', branches, 'for ', phase, 'at', gcarc)

            nbranches = len(branches)
            gc_diff = np.diff(path[:,0])
            branch_breaks = np.where(gc_diff < 0)

            #create paths if they don't exist
            try:
                grp = ds['paths'][phase][source_depth]
            except KeyError:
                ds['paths'][phase].create_group(source_depth)
            try:
                grp = ds['paths'][phase][source_depth][gcarc]
            except KeyError:
                ds['paths'][phase][source_depth].create_group(gcarc)

            #add paths for each branch
            if nbranches > 1:
                for i in range(0,nbranches):
                    branch = branches[i]

                    if i == 0:
                        bb_0 = 0
                        bb_1 = branch_breaks[0][i-1] + 1
                    else:
                        bb_0 = branch_breaks[0][i-1] + 1

                        try:
                            bb_1 = branch_breaks[0][i] + 1
                        except IndexError:
                            bb_1 = None

                    if verbose:
                        print(phase,branch,gcarc,bb_0, bb_1)

                    ds['paths'][phase][source_depth][gcarc].create_group(branch)
                    ds['paths'][phase][source_depth][gcarc][branch].create_dataset('radius', data=path[:,3][bb_0:bb_1])
                    ds['paths'][phase][source_depth][gcarc][branch].create_dataset('distance', data=path[:,4][bb_0:bb_1])
                    ds['paths'][phase][source_depth][gcarc][branch].create_dataset('time', data=path[:,7][bb_0:bb_1])

            elif nbranches == 1:

                ds['paths'][phase][source_depth][gcarc].create_group(branches[0])
                ds['paths'][phase][source_depth][gcarc][branches[0]].create_dataset('radius', data=path[:,3])
                ds['paths'][phase][source_depth][gcarc][branches[0]].create_dataset('distance', data=path[:,4])
                ds['paths'][phase][source_depth][gcarc][branches[0]].create_dataset('time', data=path[:,7])

            else:
                raise ValueError('something went wrong finding branches')

def get_travel_times1D(table,distance_in_degree,source_depth_in_km,phase,
                       component='PSV',branch=None):

    '''
    Get 1D travel travel times from a lookup table

    params:
    table <str>: path to hdf5 travel time table
    distance_in_degree <float>: source/receiver distance in degrees
    source_depth_in_km <float>: source depth in km
    phase <str>: seismic phase
    component <str>: 'PSV' or 'SH'. defaults to 'PSV'(for now)
    branch <str>: branch (doesn't do anything yet)

    returns
    time_in_s <float>: travel time in seconds.
    '''

    if branch == None and phase != 'PKP':
        branch = '1'  #defaults to branch 1 if not given
    elif branch == None and phase == 'PKP':
        branch = 'ab' #defaults to ab branch if not given

    branch = branch.encode('utf-8')
    ttt = h5py.File(table)

    #perform checks---------------------------------------------------------
    if phase not in ttt['travel_times']['1D'][component].keys():
        raise ValueError('phase {} not present in the {} table'.format(phase,table))

    #TODO add attributes to the table such as min/max depth and distance
    evdp_max = 670.0 #this is just a stand-in value
    if source_depth_in_km > evdp_max:
        raise ValueError('source depth {}-km exceeds limit of {}-km'.format(source_depth_in_km,evdp_max))

    #get closest depths
    dz_table = 1
    z1 = np.int(source_depth_in_km)
    z2 = z1 + dz_table

    #get datasets at each source depth
    dset1 = ttt['travel_times']['1D'][component][phase][str(z1)]
    dset2 = ttt['travel_times']['1D'][component][phase][str(z2)]
    b_inds1 = np.where(dset1['branch'].value==branch)[0]
    b_inds2 = np.where(dset2['branch'].value==branch)[0]
    b_inds1 = list(b_inds1)
    b_inds2 = list(b_inds2)

    if len(b_inds1) == 0 or len(b_inds2) == 0:
        raise ValueError('distance {} for phase {} and branch {} not found'.format(distance_in_degree,phase,branch))

    #get interpolated values for each depth
    ttz1 = np.interp(distance_in_degree,
                     #dset1['distance_in_degree'][:],
                     dset1['distance_in_degree'][b_inds1],
                     #dset1['time'][:],
                     dset1['time'][b_inds1])
    ttz2 = np.interp(distance_in_degree,
                     #dset2['distance_in_degree'][:],
                     dset2['distance_in_degree'][b_inds2],
                     #dset2['time'][:],
                     dset2['time'][b_inds2])

    #interp between travel times for different source depths
    w1 = (z2 - source_depth_in_km) / dz_table
    w2 = (source_depth_in_km - z1) / dz_table
    time_in_s = ttz1 * w1 + ttz2 * w2

    return time_in_s

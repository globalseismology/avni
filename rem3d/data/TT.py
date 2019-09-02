import os
import h5py
import glob
import numpy as np
import pandas as pd
from itertools import islice

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

    source_depth = str(source_depth)

    #open hdf5 file
    if os.path.exists(tt_table_name):
        ds = h5py.File(tt_table_name,'r+')
    else: 
        ds = h5py.File(tt_table_name,'w')

        #create h5py groups for travel times and paths
        ds.create_group('travel_times')
        ds['travel_times'].create_group('1D')
        ds['travel_times']['1D'].create_group(component)
        ds.create_group('paths')

    
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
            path = np.loadtxt(path_) 

            #get branch info
            with open(path_) as file_:
                branches = []
                for line_ in file_:
                    if line_.startswith("#BRANCH"):
                        branches.append(line_.strip().split()[1])

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
                    ds['paths'][phase][source_depth][gcarc][branch].create_dataset('time', data=path[:,5][bb_0:bb_1])

            elif nbranches == 1:

                ds['paths'][phase][source_depth][gcarc].create_group(branches[0])
                ds['paths'][phase][source_depth][gcarc][branches[0]].create_dataset('radius', data=path[:,3])
                ds['paths'][phase][source_depth][gcarc][branches[0]].create_dataset('distance', data=path[:,4])
                ds['paths'][phase][source_depth][gcarc][branches[0]].create_dataset('time', data=path[:,5])

            else:
                raise ValueError('something went wrong finding branches')

def get_travel_times1D(distance_in_degree,source_depth_in_km,phase,
                       component='PSV',model='NREM1D',branch='1'):

    '''
    Get 1D travel travel times from a lookup table

    params:
    distance_in_degree <float>: source/receiver distance in degrees
    source_depth_in_km <float>: source depth in km
    phase <str>: seismic phase
    component <str>: 'PSV' or 'SH'. defaults to 'PSV'(for now)
    model <str>: name of model. defaults to NREM1D (only option atm)
    branch <str>: branch (doesn't do anything yet)

    returns
    time_in_s <float>: travel time in seconds. 
    '''

    #NOTE the models below are incomplete. Once they are complete I will move
    #     them to a shared location
    if model=='NREM1D':
        ttt = h5py.File('/home/rmaguire/MODELS/NREM1D_draft_tt_table.h5','r')
    elif model=='PREM':
        ttt = h5py.File('/home/rmaguire/MODELS/PREM_tt_table.h5','r')

    #perform checks---------------------------------------------------------
    if phase not in ttt['travel_times']['1D'][component].keys():
        raise ValueError('phase {} not present in the {} table'.format(phase,model))

    #TODO add attributes to the table such as min/max depth and distance
    evdp_max = 100.0 #this is just a stand-in value 
    if source_depth_in_km > evdp_max:
        raise ValueError('depth {} is outside of available range'.format(evdp_max))

    #get closest depths
    dz_table = 1
    z1 = np.int(source_depth_in_km)
    z2 = z1 + dz_table

    #get datasets at each source depth
    dset1 = ttt['travel_times']['1D'][component][phase][str(z1)]
    dset2 = ttt['travel_times']['1D'][component][phase][str(z2)]

    #get interpolated values for each depth
    ttz1 = np.interp(distance_in_degree,
                     dset1['distance_in_degree'][:],
                     dset1['time'][:])
    ttz2 = np.interp(distance_in_degree,
                     dset2['distance_in_degree'][:],
                     dset2['time'][:])

    print(ttz1,ttz1.shape)
    print(ttz2,ttz2.shape)

    #interp between travel times for different source depths
    w1 = (source_depth_in_km - z1) / dz_table
    w2 = (z2 - source_depth_in_km) / dz_table
    time_in_s = ttz1 * w1 + ttz2 * w2 

    return time_in_s

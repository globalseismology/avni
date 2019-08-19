import os
import h5py
import glob
import numpy as np
from itertools import islice

def write_h5_tt_table(cagc_rays_outdir,tt_table_name,source_depth,verbose=True):
    '''
    writes an hdf5 file based on the output of CAGCRAYS

    params:
    cagc_rays <str>: path to output directory of CAGCRAYS
    tt_table_name <str>: name of hdf5 file to create or add to
    '''

    source_depth = str(source_depth)

    #open hdf5 file
    if os.path.exists(tt_table_name):
        ds = h5py.File(tt_table_name,'r+')
    else: 
        ds = h5py.File(tt_table_name,'w')
    
    #get list of directories containing travel time and path into
    phase_list = glob.glob(cagc_rays_outdir+'/*')
     
    #create h5py groups for travel times and paths
    ds.create_group('travel_times')
    ds.create_group('paths')

    for phase_ in phase_list:
        phase = phase_.split('/')[-1]

        #get travel times
        try:
            tts = np.genfromtxt(phase_+'/ttandtstar.txt',skip_header=1,dtype=None)
            tts = np.flipud(tts) #arange in order of increasing distance
        except(OSError):
            raise ValueError('The directory ',phase,'does not contain ttandtstar.txt')

        #create a new group for each phase
        try:
            ds['travel_times'].create_group(phase) 
            ds['travel_times'][phase].create_group(source_depth)
        except RuntimeError:
            print('group already exists... appending new data to it')

        #add data for each phase (when applicable follow obspy.TauP naming)
        ds['travel_times'][phase][source_depth].create_dataset('distance_in_degree',data=tts['f0'])
        ds['travel_times'][phase][source_depth].create_dataset('ray_param_sec_degree',data=tts['f1'])
        ds['travel_times'][phase][source_depth].create_dataset('time',data=tts['f2'])
        ds['travel_times'][phase][source_depth].create_dataset('dDdp',data=tts['f3'])
        ds['travel_times'][phase][source_depth].create_dataset('tstar',data=tts['f4'])
        ds['travel_times'][phase][source_depth].create_dataset('branch',data=tts['f5'])

        #get raypaths
        path_list = glob.glob(cagc_rays_outdir+'/'+phase+'/path_SYNTH*')
        path_list.sort()

        try:
            ds['paths'].create_group(phase)
        except RuntimeError:
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

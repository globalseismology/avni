#!/usr/bin/env python

#######################################################################################
# python 3 compatibility
from __future__ import absolute_import, division, print_function

#####################  IMPORT STANDARD MODULES   ######################################

import sys,os
import six
import numpy as np #for numerical analysis
import fortranformat as ff #reading/writing fortran formatted text
import pandas as pd
import h5py
import time
import warnings
from .. import constants

if sys.version_info[0] >= 3: unicode = str

####################       I/O ROUTINES     ######################################

def get_travel_times1D(table,distance_in_degree,period,output='pvel',overtone=0,phase=None):
    """
    Return the arrival time in seconds of a surface-wave phase or mode type/overtone.
    When phase is queries, mode_type/overtone/arc are ignored

    Input Parameters:
    ----------------

    table <str>: path to hdf5 format modes table
    distance_in_degree <float>: great circle distance, can be minor or major arc
    period <float>: period in second
    mode_type <str>: either "spheroidal", "toroidal", or "radial"
    overtone <int>: overtone number (defaults to 0, ie. fundamental mode)

    """

    # Perform some checks
    if phase == None: #if no phase is given, assume minor arc
        mode_type='spheroidal'
        iorbit=1
    else:

        if isinstance(phase, str):
            # heuristics: 'R1' - call spheroidal, G1=toroidal
            if phase.startswith('L') or phase.startswith('G'):
                mode_type='toroidal'
            elif phase.startswith('R'):
                mode_type='spheroidal'
            iorbit=int(phase[1])
        else:
            raise ValueError('Only phases like R1, G1 can be called')


    # Perform interpolation to get velocity for requested period
    vel_i = get_velocity(table=table,period=period,overtone=overtone,mode_type=mode_type,output=output)

    # Find travel time
    if iorbit==1:
        distance_in_km = constants.deg2km.magnitude * distance_in_degree
        igc_count=0
    else:
        if np.mod(iorbit,2) == 0: #even orbit hence G2,G4, etc.
            igc_count=iorbit/2   #number of times almost circled the earth
            distance_in_km = constants.deg2km.magnitude * (360.*float(igc_count)-distance_in_degree)
        else:
            igc_count=(iorbit-1)/2
            distance_in_km = constants.deg2km.magnitude * (360.*float(igc_count)+distance_in_degree)

    time_in_s =  distance_in_km * (1./vel_i)

    return time_in_s

def get_velocity(table,period,overtone,mode_type,output='pvel'):
    """
    Return the arrival time in seconds of a surface-wave phase or mode type/overtone.
    When phase is queries, mode_type/overtone/arc are ignored

    Input Parameters:
    ----------------

    table <str>: path to hdf5 format modes table
    distance_in_degree <float>: great circle distance, can be minor or major arc
    period <float>: period in second
    mode_type <str>: either "spheroidal", "toroidal", or "radial"
    overtone <int>: overtone number (defaults to 0, ie. fundamental mode)

    """

    table = h5py.File(table,'r')
    omega_query = (1./period) * (2.*np.pi) #the requested frequency in rad/s

    if mode_type=='S':
        mode_type='spheroidal'
    elif mode_type=='T':
        mode_type='toroidal'
    elif mode_type=='R':
        mode_type='radial'

    omega = table[mode_type][str(overtone)].attrs['omega']
    vel = table[mode_type][str(overtone)].attrs[output]

    if omega_query < np.min(omega) or omega_query > np.max(omega):
        raise ValueError('period {} doesnt exist in table for requested mode'.format(period))

    # Perform interpolation to get velocity for requested period
    vel_i = np.interp(omega_query,omega,vel)

    return vel_i

def get_dispersion_curve(table,mode_type,output='pvel',overtone=0,freq_units='mhz'):
    '''
    return a dispersion curve for

    params:
    table <str>: path to hdf5 format modes table
    mode_type <str>: either "spheroidal", "toroidal", or "radial"
    output <str>: either "gvel" or "pvel" for group or phase velocity
    overtone <int>: overtone number (defaults to 0, ie. fundamental mode)
    freq_units <str>: units of frequency axis. either "rad/s","hz",or "mhz"

    returns:
    freq,vel: frequency (in freq_units) and (group or phase) velocity in km/s
    '''

    table = h5py.File(table,'r')
    omega = table[mode_type][str(overtone)].attrs['omega']
    vel = table[mode_type][str(overtone)].attrs[output]

    if freq_units.lower() == 'hz':
        freq = omega/(2.*np.pi)
    elif freq_units.lower() == 'mhz':
        freq = (omega/(2.*np.pi)) * 1000.
    elif freq_units.lower() == 'rad/s':
        freq = omega

    return freq,vel

def readSWascii(file, delim = '-',required = None,warning=False):
    """Reads the AVNI format for analysis and plotting.

    Input parameters:
    ----------------

    file :  input file in default AVNI format

    delim : delimiter that combines fields into a joint field e.g. network-station
            seperate out during I/O.

    required : fields needed as comments (e.g. #CITE: Author et al., YEAR) in the file
                defaults: 'CITE', 'SHORTCITE', 'REFERENCE MODEL', 'PVEL', 'CRUST',
                          'MODEL3D', 'SIGMATYPE', 'WEITYPE','EQTYPE', 'STATTYPE',
                          'FORMAT', 'WRITE', 'FIELDS'

    Output:
    ------

    SWdata :  dictionary with fields data, metadata and comments
    """
    # defaults
    if required == None: required = ['CITE', 'SHORTCITE', 'REFERENCE MODEL', 'PVEL', 'CRUST', 'MODEL3D', 'SIGMATYPE', 'WEITYPE','EQTYPE', 'STATTYPE', 'FORMAT', 'WRITE', 'FIELDS']
    if (not os.path.isfile(file)): raise IOError("Filename ("+file+") does not exist")

    # checks for CITE and SHORTCITE comments
    comments = []
    for line in open(file):
        if line.startswith("#"): comments.append(line.rstrip())

    metadata={}; notes=[]
    for line in comments:
        if line.startswith("#") and ':' in line and not line.startswith("#NOTES:"):
            key = line[1:].split(':')[0]
            value = line[1:].split(':')[1].rstrip().lstrip()
            metadata[key] = value
        elif line.startswith("#NOTES:"):
            notes.append(line)
    metadata['COMMENTS'] = notes

    #check if required fields exist
    for key in required:
        try:
            namelist = metadata[key].split()
        except KeyError:
            if warning:
                print(key+" should be defined as a comment in "+file)
            else:
                raise KeyError(key+" should be defined as a comment in "+file)

    data = pd.read_table(file,header=None,comment='#',sep='\s+',names=metadata['FIELDS'].split())

    # replace the fields by divided fields if concatenated
    for column in data.columns.tolist():
        if delim in column:
            newcolumns = column.split(delim)
            try:
                data[newcolumns] = data[column].str.split(delim,expand=True)
                data = data.drop(column, 1) # get rid of old composite column
                data.replace('', np.nan, inplace=True) #replace blank with nan
            except:
                warnings.warn('could not split the column '+column+' with delimiter '+delim)
    if 'stat' in data.columns.values and 'net' in data.columns.values:
        data['path'] = data['cmtname'] + '_'+ data['stat'] + '-' + data['net']

    SWdata = {}; SWdata['data'] = data; SWdata['metadata'] = metadata
    return SWdata

def writeSWascii(SWdata,filename,iflagthreshold=None,delim='-',writeheader=True,writedata=True,verbose=True):
    """
    Writes the AVNI format for analysis and plotting

    Input parameters:
    ----------------

    SWdata : dictionary of SW data with metadata and data fields

    filename : output file name

    delim : delimiter that combines fields into a joint field e.g. network-station
            seperate out during I/O.

    iflagthreshold : threshold for iflag which corresponds to the processing level
                     that was cleared

    """


    metadata = SWdata['metadata']
    header_line  =  ff.FortranRecordWriter('('+metadata['WRITE']+')')

    # write out header
    if writeheader:
        printstr  =  [unicode('#'+key+':'+metadata[key]+'\n') for key in metadata.keys() if key not in ['FIELDS','WRITE','COMMENTS']]
        printstr.append(unicode('#'*len(metadata['FIELDS'])+'\n'))
        for comment in metadata['COMMENTS']: printstr.append(unicode(comment+'\n'))
        printstr.append(unicode('#'*len(metadata['FIELDS'])+'\n'))
        printstr.append(unicode('#WRITE:'+metadata['WRITE']+'\n'))
        printstr.append(unicode('#FIELDS:'+metadata['FIELDS']+'\n'))

        with open(filename,'w') as f: f.writelines(printstr)

    # append the series to output file (one column per row now)
    if writedata:
        if iflagthreshold is None:
            data = SWdata['data']
        else:
            data = SWdata['data'][SWdata['data']['iflag'] >= iflagthreshold]

        if len(data) == 0: raise ValueError('No data row clears threshold for ', iflagthreshold)
        if isinstance(data, pd.Series): data = pd.DataFrame(data.to_dict(), index=[0])
        namelist = metadata['FIELDS'].split()
        # replace the fields by divided fields if concatenated
        for column in namelist:
            if delim in column and column not in data.columns.values:
                oldcolumns = column.split(delim)

                data[column]=data[oldcolumns].fillna('').apply(lambda x: delim.join(x.astype(str).values), axis=1)
                data = data.drop(oldcolumns, 1) # get rid of old separate columns
        # reorder according to FIELDS
        data = data.reindex(namelist,axis=1)

        # all pandas version
        cols=metadata['FIELDS']

        # initialize fortranformat writer -- apply by row
        line = ff.FortranRecordWriter('('+metadata['WRITE']+')')

        # using apply
        for field in ['cmtname','stat-net-chan-loc']:
            data[field] = data[field].str.ljust(15)
        Formated_Series=data.apply(lambda x : line.write(x.values),axis=1)
        Formated_Series.to_csv(filename,index=False,header=False,mode='a')
    if verbose: print("....written "+str(len(data))+" observations to "+filename)
    return

def SWasciitohdf5(files,hdffile = 'Summary.SW.data.h5',datatype='summary',delim='-'):
    """
    Read a list of files to hdf5 container format.

    Input parameters:
    ----------------

    files: a panda list of ascii files to append

    hdffile: output hdf5 file to append fields to

    datatype: type of data to group the data with

    delim : delimiter that combines fields into a joint field e.g. network-station
            seperate out during I/O.

    """
    hf = h5py.File(hdffile, 'a')
    for _ , row in files.iterrows():
        file = row[0]
        print("... creating group folders based on "+file+" in "+hdffile+" with the "+datatype+" group")
        SWdata = readSWascii(file,delim=delim)
        metadata = SWdata['metadata']; data = SWdata['data']
        group = metadata['SHORTCITE']
        uniquefolders = ['overtone', 'period']
        for ii in np.arange(len(uniquefolders)):
            key = uniquefolders[ii]
            if len(np.unique(data[key])) != 1: raise ValueError('Number  of unique values of '+key+' should be 1 in '+file)
            if ii == 0:
                folder = str(np.unique(data[key])[0])
            else:
                folder = folder+'/'+str(np.unique(data[key])[0])
            g1 = hf.require_group(folder)
            g1.attrs['name']=key
            if key == 'period': g1.attrs['units']='s'

        # Loop over every type or orbit and store values e.g. R1,R2,R3
        for typeiorb in np.unique(data['typeiorb']):
            folderorb = folder+'/'+typeiorb
            g1 = hf.require_group(folderorb)
            g1.attrs['name']='typeiorb'
            g2 = g1.require_group(group)
            g2.attrs['FILE'] = file
            for field in metadata.keys():
                if field == 'COMMENTS':
                    asciiList = [n.encode("ascii", "ignore") for n in metadata['COMMENTS']]
                    g2.attrs.create('COMMENTS', asciiList, (len(asciiList),1),'S200')
                else:
                    g2.attrs[field]=metadata[field]
            g3 = g2.require_group('data/'+datatype)
            subfolder = folderorb+'/'+group+'/data/'+datatype

            #### store path ####
            g4 = g3.require_group('arrival/paths')
            g4.attrs['desc'] = 'ipath is the unique path in CMT_STATION-NETWORK format'

            #### store source ####
            g4 = g3.require_group('source/'+metadata['EQTYPE'])

            #### store station  ####
            g4 = g3.require_group('station/'+metadata['STATTYPE'])

            #### store data and predictions  ####
            g4 = g3.require_group('arrival/observed')
            g4.attrs['units']='s'
            g4.attrs['desc'] = 'observed relative arrivals w.r.t. various 1D reference models'
            g5 = g4.require_group(metadata['REFERENCE MODEL'])
            g5.attrs['pvel'] = float(metadata['PVEL'].split()[0])
            g5.attrs['pvel_units'] = metadata['PVEL'].split()[1]

            g4 = g3.require_group('arrival/predicted')
            g4.attrs['units']='s'
            g4.attrs['desc'] = 'predicted relative and absolute arrivals based on various 1D and 3D mantle and crustal models'

            # store weights/uncertainties
            g4 = g3.require_group('arrival/weight')
            g4.attrs['desc'] = 'weight of the measurement to account fo uneven sampling etc.'
            g4 = g3.require_group('arrival/sigma')
            g4.attrs['desc'] = 'uncertainity of the measurement'

            #### store flag ####
            g4 = g3.require_group('arrival/others')
            g4.attrs['desc'] = 'other fields relevant to the arrival e.g. iflag describes the processing scheme : 0 is raw'

    hf.close()

    ####################################################
    # Now store relevant pandas Dataframe

    for _ , row in files.iterrows():
        file = row[0]
        print("... adding "+file+" to "+hdffile+" in the "+datatype+" group")
        SWdata = readSWascii(file,delim=delim)
        metadata = SWdata['metadata']; data = SWdata['data']
        group = metadata['SHORTCITE']
        uniquefolders = ['overtone', 'period']
        for ii in np.arange(len(uniquefolders)):
            key = uniquefolders[ii]
            if len(np.unique(data[key])) != 1: raise ValueError('Number  of unique values of '+key+' should be 1 in '+file)
            if ii == 0:
                folder = str(np.unique(data[key])[0])
            else:
                folder = folder+'/'+str(np.unique(data[key])[0])

        # Loop over every type or orbit and store values e.g. R1,R2,R3
        for typeiorb in np.unique(data['typeiorb']):
            folderorb = folder+'/'+typeiorb
            subfolder = folderorb+'/'+group+'/data/'+datatype

            #### store path ####
            columns = ['distkm', 'path','delellphase']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/arrival/paths'),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store source ####
            columns = ['cmtname', 'eplat', 'eplon', 'epdep']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/source/'+metadata['EQTYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store station  ####
            columns = ['stlat', 'stlon','stat','net', 'chan', 'loc']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/station/'+metadata['STATTYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store data and predictions  ####
            columns = ['refphase', 'delobsphase']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/arrival/observed/'+metadata['REFERENCE MODEL']),format='table', mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            if metadata['MODEL3D'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delpredphase'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/predicted/'+metadata['REFERENCE MODEL']+'/'+metadata['MODEL3D']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
            if metadata['CRUST'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delcruphase'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/predicted/'+metadata['REFERENCE MODEL']+'/'+metadata['CRUST']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            # store weights/uncertainties
            if metadata['WEITYPE'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delweight'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/weight/'+metadata['WEITYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
            if metadata['SIGMATYPE'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delsigma'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/sigma/'+metadata['SIGMATYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store flag ####
            try:
                data.loc[data['typeiorb'] == typeiorb]['iflag'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/others'),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
            except KeyError:
                if datatype == 'raw':
                    iflag = pd.DataFrame(0, index=np.array(data.loc[data['typeiorb'] == typeiorb].index), columns=['iflag'])
                    iflag.to_hdf(hdffile, key=unicode(subfolder+'/arrival/others'),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
                else:
                    raise ValueError('iflag not defined for datatype ('+datatype+') unless specified for each measurement row in '+file)


def SWhdf5toascii(query = '0/25.0/L1/REM3D',hdffile = 'Summary.SW.data.h5',iflag=0,datatype='summary',delim='-', outfile=None,model3d = None, refmodel = None,crust=None,weitype=None, sigmatype= None,stattype =None,eqtype=None):
    """
    write hdf field to a file. None is the default i.e. the values in metadata

    Input parameters:
    ----------------

    query: key query to the hdf5 file till the group level e.g. 0/25.0/L1/GDM52

    hdffile: output hdf5 file to append fields to

    datatype: type of data to group the data with e.g. raw, processed, summary

    delim : delimiter that combines fields into a joint field e.g. network-station
            seperate out during I/O.

    iflag : select the flag to sub-select data

    outfile : output ascii file. If None, decided based on query

    model3d, refmodel,crust,weitype :   Fields used to query various groups and attributes
    sigmatype,stattype ,eqtype      :   If None, use the values in the hdffile metadata

    Output:
    ------

    ASCII file with the values selected

    """
    if (not os.path.isfile(hdffile)): raise IOError("Filename ("+hdffile+") does not exist")
    # read the dataframes from hdf5
    SWdata = readSWhdf5(query,hdffile,iflag,datatype,delim,model3d,refmodel,crust, weitype, sigmatype,stattype,eqtype)

    # write the pandas to file
    fields = query.split('/')
    if outfile is None: filename = fields[3]+'_'+fields[0]+'_'+fields[1]+'_'+fields[2]+'.AVNI'
    writeSWascii(SWdata,filename,delim=delim)

def readSWhdf5(query = '0/25.0/L1/REM3D',hdffile = 'Summary.SW.data.h5',iflag=0,datatype='summary',delim='-',model3d=None, refmodel=None,crust=None,weitype=None, sigmatype= None,stattype =None,eqtype=None):
    """
    Read a list of files to hdf5 container format.

    Input parameters:
    ----------------

    query: key query to the hdf5 file till the group level e.g. 0/25.0/L1/GDM52

    hdffile: output hdf5 file to append fields to

    datatype: type of data to group the data with e.g. raw, processed, summary

    delim : delimiter that combines fields into a joint field e.g. network-station
            seperate out during I/O.

    iflag : select the flag to sub-select data

    model3d, refmodel,crust,weitype :   Fields used to query various groups and attributes
    sigmatype,stattype ,eqtype      :   If None, use the values in the hdffile metadata

    Output:
    ------

    SWdata : dictionary of SW data with metadata and data fields consistent with
             readSWascii and writeSWascii
    """
    SWdata = {}
    if (not os.path.isfile(hdffile)): raise IOError("Filename ("+hdffile+") does not exist")
    hf = h5py.File(hdffile, 'r')
    SWdata['metadata'] = {}
    for name,value in hf[query].attrs.items():
        if name == 'COMMENTS':
            SWdata['metadata'][name] = [txt[0].decode('utf-8') for txt in hf[query].attrs['COMMENTS'][:].tolist()]
        else:
            SWdata['metadata'][name]=value
    try:
        #### store flag ####
        df0 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/arrival/others',where=['iflag=='+str(iflag)])
        #### get path ####
        df1 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/arrival/paths', where=df0.index)
        #### get source ####
        if eqtype is None: eqtype = SWdata['metadata']['EQTYPE']
        df2 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/source/' +eqtype, where=df0.index)
        #### store station  ####
        if stattype is None: stattype = SWdata['metadata']['STATTYPE']
        df3 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/station/' +stattype, where=df0.index)
        #### store data and predictions  ####
        if refmodel is None: refmodel = SWdata['metadata']['REFERENCE MODEL']
        df4 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/arrival/observed/' +refmodel, where=df0.index)
        if model3d is None: model3d = SWdata['metadata']['MODEL3D']
        if model3d == 'None':
            df5 = pd.DataFrame(0, index=np.array(df3.index), columns=['delpredphase'])
        else:
            df5 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/arrival/predicted/' +refmodel+'/'+model3d, where=df0.index)
        if crust is None: crust = SWdata['metadata']['CRUST']
        if crust == 'None':
            df6 = pd.DataFrame(0, index=np.array(df3.index), columns=['delcruphase'])
        else:
            df6 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/arrival/predicted/' +refmodel+'/'+crust, where=df0.index)

        # store weights/uncertainties
        if weitype is None: weitype = SWdata['metadata']['WEITYPE']
        if weitype == 'None':
            df7 = pd.DataFrame(0, index=np.array(df3.index), columns=['delweight'])
        else:
            df7 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/arrival/weight/' +weitype, where=df0.index)
        if sigmatype is None: sigmatype = SWdata['metadata']['SIGMATYPE']
        if sigmatype == 'None':
            df8 = pd.DataFrame(0, index=np.array(df3.index), columns=['delsigma'])
        else:
            df8 = pd.read_hdf(hdffile, query+'/data/'+datatype+'/arrival/sigma/' +sigmatype, where=df0.index)

        SWdata['data'] = pd.concat([df0,df1,df2,df3,df4,df5,df6,df7,df8],axis=1)
    except:
        raise ValueError('fields missing in location '+query+'/data/'+datatype)

    # write the pandas to file
    fields = query.split('/')
    SWdata['data']['overtone'] = int(fields[0])
    SWdata['data']['period'] = float(fields[1])
    SWdata['data']['typeiorb'] = fields[2]

    # update the metadata attributes based on input query
    SWdata['metadata']['REFERENCE MODEL'] = refmodel
    findindx = query+'/data/'+datatype+'/arrival/observed/' +refmodel
    SWdata['metadata']['PVEL'] = str(hf[findindx].attrs.__getitem__('pvel'))+' '+hf[findindx].attrs.__getitem__('pvel_units')
    SWdata['metadata']['REFERENCE MODEL'] = refmodel
    SWdata['metadata']['CRUST'] = crust
    SWdata['metadata']['MODEL3D'] = model3d
    SWdata['metadata']['SIGMATYPE'] = sigmatype
    SWdata['metadata']['WEITYPE'] = weitype
    SWdata['metadata']['EQTYPE'] = eqtype
    SWdata['metadata']['STATTYPE'] = stattype
    return SWdata

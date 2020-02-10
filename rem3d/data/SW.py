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
import progressbar
import pdb

if sys.version_info[0] >= 3: unicode = str

####################       I/O ROUTINES     ######################################

def readSWascii(file, delim = '-',required = None,warning=False):
    """Reads the REM3D format for analysis and plotting.

    Input parameters:
    ----------------

    file :  input file in default REM3D format

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
            data[newcolumns] = data[column].str.split(delim,expand=True)
            data = data.drop(column, 1) # get rid of old composite column
            data.replace('', np.nan, inplace=True) #replace blank with nan
    data['path'] = data['cmtname'] + '_'+ data['stat'] + '-' + data['net']

    SWdata = {}; SWdata['data'] = data; SWdata['metadata'] = metadata
    return SWdata

def writeSWascii(SWdata,filename,iflagthreshold=None,delim='-'):
    """
    Writes the REM3D format for analysis and plotting

    Input parameters:
    ----------------

    SWdata : dictionary of SW data with metadata and data fields

    filename : output file name

    delim : delimiter that combines fields into a joint field e.g. network-station
            seperate out during I/O.

    iflagthreshold : threshold for iflag which corresponds to the processing level
                     that was cleared

    """

    if iflagthreshold is None:
        data = SWdata['data']
    else:
        data = SWdata['data'][SWdata['data']['iflag'] >= iflagthreshold]

    metadata = SWdata['metadata'];
    header_line  =  ff.FortranRecordWriter('('+metadata['WRITE']+')')

    printstr  =  [unicode('#'+key+':'+metadata[key]+'\n') for key in metadata.keys() if key not in ['FIELDS','WRITE','COMMENTS']]
    printstr.append(unicode('#'*len(metadata['FIELDS'])+'\n'))
    for comment in metadata['COMMENTS']: printstr.append(unicode(comment+'\n'))
    printstr.append(unicode('#'*len(metadata['FIELDS'])+'\n'))
    for key in metadata.keys():
        if key in ['FIELDS','WRITE']: printstr.append(unicode('#'+key+':'+metadata[key]+'\n'))

    namelist = metadata['FIELDS'].split()
    # replace the fields by divided fields if concatenated
    for column in namelist:
        if delim in column:
            oldcolumns = column.split(delim)

            data[column]=data[oldcolumns].fillna('').apply(lambda x: delim.join(x.astype(str).values), axis=1)
            data = data.drop(oldcolumns, 1) # get rid of old separate columns
    # reorder according to FIELDS
    data = data.reindex(namelist,axis=1)

    f = open(filename,'w')
    f.writelines(printstr)
    for ii in progressbar.progressbar(range(len(data))):
        line=[]
        for val in data.values[ii]:
            if isinstance(val,six.string_types):
                line.append(val.ljust(15))
            else:
                line.append(val)
        try:
            arow = header_line.write(line)
        except:
            raise IOError('No line to print')
        f.write(unicode(arow+'\n'))
    f.close()
    print("....written "+str(len(data))+" observations to "+filename)
    return


def SWasciitohdf5(files,hdffile = 'SW.rem3d.data.h5',datatype='raw',delim='-'):
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
        print("... adding "+file+" to "+hdffile+" in the "+datatype+" group")
        SWdata = readSWascii(file,delim=delim)
        metadata = SWdata['metadata']; data = SWdata['data']
        group = metadata['SHORTCITE']
        uniquefolders = ['overtone', 'period']
        for ii in np.arange(len(uniquefolders)):
            key = uniquefolders[ii]
            if len(np.unique(data[key])) is not 1: raise ValueError('Number  of unique values of '+key+' should be 1 in '+file)
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
            columns = ['distkm', 'path','delellphase']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/arrival/paths'),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store source ####
            g4 = g3.require_group('source/'+metadata['EQTYPE'])
            columns = ['cmtname', 'eplat', 'eplon', 'epdep']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/source/'+metadata['EQTYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store station  ####
            g4 = g3.require_group('station/'+metadata['STATTYPE'])
            columns = ['stlat', 'stlon','stat','net', 'chan', 'loc']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/station/'+metadata['STATTYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store data and predictions  ####
            g4 = g3.require_group('arrival/observed')
            g4.attrs['units']='s'
            g4.attrs['desc'] = 'observed relative arrivals w.r.t. various 1D reference models'
            g5 = g4.require_group(metadata['REFERENCE MODEL'])
            g5.attrs['pvel'] = float(metadata['PVEL'].split()[0])
            g5.attrs['pvel_units'] = metadata['PVEL'].split()[1]
            columns = ['refphase', 'delobsphase']
            data.loc[data['typeiorb'] == typeiorb][columns].to_hdf(hdffile, key=unicode(subfolder+'/arrival/observed/'+metadata['REFERENCE MODEL']),format='table', mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            g4 = g3.require_group('arrival/predicted')
            g4.attrs['units']='s'
            g4.attrs['desc'] = 'predicted relative and absolute arrivals based on various 1D and 3D mantle and crustal models'
            if metadata['MODEL3D'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delpredphase'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/predicted/'+metadata['REFERENCE MODEL']+'/'+metadata['MODEL3D']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
            if metadata['CRUST'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delcruphase'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/predicted/'+metadata['REFERENCE MODEL']+'/'+metadata['CRUST']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            # store weights/uncertainties
            g4 = g3.require_group('arrival/weight')
            g4.attrs['desc'] = 'weight of the measurement to account fo uneven sampling etc.'
            if metadata['WEITYPE'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delweight'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/weight/'+metadata['WEITYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
            g4 = g3.require_group('arrival/sigma')
            g4.attrs['desc'] = 'uncertainity of the measurement'
            if metadata['SIGMATYPE'] != 'None':
                data.loc[data['typeiorb'] == typeiorb]['delsigma'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/sigma/'+metadata['SIGMATYPE']),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)

            #### store flag ####
            g4 = g3.require_group('arrival/others')
            g4.attrs['desc'] = 'other fields relevant to the arrival e.g. iflag describes the processing scheme : 0 is raw'
            try:
                data.loc[data['typeiorb'] == typeiorb]['iflag'].to_hdf(hdffile, key=unicode(subfolder+'/arrival/others'),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
            except KeyError:
                if datatype == 'raw':
                    iflag = pd.DataFrame(0, index=np.array(data.loc[data['typeiorb'] == typeiorb].index), columns=['iflag'])
                    iflag.to_hdf(hdffile, key=unicode(subfolder+'/arrival/others'),format='table',mode='a',complib='bzip2',complevel=9,dropna=True,data_columns=True,append=True)
                else:
                    raise ValueError('iflag not defined for datatype ('+datatype+') unless specified for each measurement row in '+file)

    hf.close()

def SWhdf5toascii(query = '0/25.0/L1/GDM52',hdffile = 'SW.rem3d.data.h5',iflag=0,datatype='processed',delim='-', outfile=None,model3d = None, refmodel = None,crust=None,weitype=None, sigmatype= None,stattype =None,eqtype=None):
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
    if outfile is None: filename = fields[3]+'_'+fields[0]+'_'+fields[1]+'_'+fields[2]+'.REM3D'
    writeSWascii(SWdata,filename,delim=delim)

def readSWhdf5(query = '0/25.0/L1/GDM52',hdffile = 'SW.rem3d.data.h5',iflag=0,datatype='raw',delim='-',model3d=None, refmodel=None,crust=None,weitype=None, sigmatype= None,stattype =None,eqtype=None):
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

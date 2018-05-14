''' wrapper for calling rem3d commands '''
import rem3d
import numpy as np
import pandas as pd
import numpy.lib.recfunctions as nprec
import os
import processBokehData as pbd

class data_fetch(object):

    def __init__(self,AppletType):
        self.filesdir=rem3d.tools.get_filedir()
        self.AppletType=AppletType
        self.MetaData=dict()
        self.DataInitialized=False
        self.DegSigRounding=5 # for rounding station and epicenter lat/lons
        self.getDataSetList()
        print ("rem3d_wrapper initialized with AppletType " + AppletType)
        print ("rem3d file storage path is "+self.filesdir)
        return

    def pullData(self,file_name,data_id,check_for_update=True,download_if_missing=True):
        ''' will pull rem3d data into general data dictionary object '''
        if check_for_update:
            self.tryUpdate(file_name)

        if download_if_missing and not os.path.isfile(self.filesdir+'/'+file_name):
            print("file not available locally, trying to fetch it")
            self.tryUpdate(file_name)

        if not self.DataInitialized:
            self.Data=dict()
            self.MetaData=dict()
            self.DataInitialized=True

        if self.AppletType=='SurfaceWaves':
            # readREM3DSWformat returns a tuple
            # SWdata,comments,reference,shortref = rem3d.data.readREM3DSWformat(maindir+'/files/'+file)
            print("pulling "+file_name+" into memory at self.Data["+str(data_id)+"]")
            self.Data[data_id],comments,reference,shortref=rem3d.data.readREM3DSWformat(self.filesdir+'/'+file_name,use_pandas=True)
            self.MetaData[data_id]={'comments':comments,'ref':reference,'shortref':shortref,'data_id':data_id,'file_name':file_name}
            print("data pulled and stored")
        return

    def processData(self):
        ''' processes data into format for plotting '''
        if self.AppletType=='SurfaceWaves':
            self.buildCommonSWStations() # finds unique station-eq paths
            self.buildMapData() # for station, eq map
            self.compareSWStudies() # for group 2 vs group 1
        return

    def compareSWStudies(self):
        self.SW_StatEpPairs,binned_vals=rem3d.tools.compareSWStudiesPandas(self.SW_StatEpPairs)
        self.Scatter3Src=pbd.setSWScatter3(self.SW_StatEpPairs)
        return

    def buildCommonSWStations(self):
        self.SW_StatEpPairs,self.SW_stats,self.SW_eps=rem3d.tools.getcommonSWcatalogsPandas(self.Data[1],self.Data[2],self.DegSigRounding)
        return

    def buildMapData(self):
        ''' builds station/epicenter/other data for mapping '''
        if self.AppletType=='SurfaceWaves':
            try:
                self.SW_StatEpPairs
            except:
                self.buildCommonSWStations()
            print ('building map data')
            self.MapData=pbd.buildSWMapData(self.SW_stats,self.SW_eps)
        return

    def tryUpdate(self,file_name):
        rem3d.data.update_file(self.filesdir,file_name)
        return

    def updateSelected(self,inds):
        if len(inds) == 0:
            Nstats=len(self.Data)
            selectednets = list(set(self.Data['data_id'])) # pulls unique networks
        else:
            Nstats=len(inds)
            selectednets = list(set(self.Data['data_id'])) # pulls unique networks

        selectednets.sort()
        stattext="Selected Station Details\n\nStation Count: "+str(Nstats)+"\nReferences Selected:"
        for net in selectednets:
            stattext=stattext+"\n    "+net
        stattext=stattext+"\n\nSlider: depth\nSlider:?\nButton: Calculate\nButton: Download"
        return stattext;

    def getDataSetList(self):
        ''' builds list of available data '''

        if self.AppletType=='SurfaceWaves':
            #GDM52.0.R1.100s.interp.rem3d --> study.mode.arcpath.period.type.extension
            self.DataSetList=pd.DataFrame(columns=['study','mode','path','period_s','freq_Hz','file_name'])
            file_list=os.listdir(self.filesdir)
            file_list.sort()
            col_list=['study','normalmode','path','period_s','freq_Hz','file_name']
            for fi in file_list:
                if len(fi.split('.'))==6:
                    per=int(fi.split('.')[3].rstrip('s'))
                    rowdata=[[fi.split('.')[0],int(fi.split('.')[1]),fi.split('.')[2],per,1.0/per,fi]]
                    self.DataSetList=self.DataSetList.append(pd.DataFrame(rowdata,columns=col_list))
            # print self.DataSetList
        return

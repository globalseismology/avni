''' wrapper for calling rem3d commands '''
import rem3d
import numpy as np
import numpy.lib.recfunctions as nprec
import os

class data_fetch(object):

    def __init__(self,AppletType):
        relapath='../rem3d/rem3d'
        self.rem3ddir=os.path.dirname(os.path.abspath(relapath))
        self.AppletType=AppletType
        self.MetaData=dict()
        self.DataInitialized=False
        print ("rem3d_wrapper initialized with AppletType " + AppletType)
        print ("rem3d path is "+self.rem3ddir)
        return

    def pullData(self,file_name,data_id,check_for_update=True):
        ''' will pull rem3d data into general data dictionary object '''
        if check_for_update:
            self.tryUpdate(file_name)

        if self.AppletType=='SurfaceWaves':
            # readREM3DSWformat returns a tuple
            # SWdata,comments,reference,shortref = rem3d.data.readREM3DSWformat(maindir+'/files/'+file)
            print("pulling "+file_name+" into memory at self.Data["+data_id+"]")
            SWdata,comments,reference,shortref=rem3d.data.readREM3DSWformat(self.rem3ddir+'/files/'+file_name)
            self.MetaData[data_id]=(comments,reference,shortref,file_name)
            # append an id flag
            SWdata=nprec.append_fields(SWdata,'data_id',np.array([data_id for _ in range(len(SWdata))]),usemask=False)
            # SWdata=np.insert(SWdata, 0, data_id, axis=1)
            if self.DataInitialized:
                self.Data=np.concatenate((self.Data,SWdata),axis=0)
            else:
                self.Data=SWdata
                self.DataInitialized=True
        return

    def buildMapData(self):        
        if self.AppletType=='SurfaceWaves':
            self.MapData=np.unique(self.Data[['stat','stlat','stlon','data_id']])
        return

    def tryUpdate(self,file_name):
        rem3d.data.update_file(file_name)
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

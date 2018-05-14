''' methods for processing data for bokeh plots '''

from bokeh.plotting import ColumnDataSource
from bokeh.colors import RGB
from bokeh.models.widgets import PreText
import pandas as pd

def buildSWMapData(Stats,Eps):

    SrcDict=dict()
    SrcDict['title']='Stations and Epicenters' # title of scatter plots
    SrcDict['nsets']=2 # number of scatter plots

    # put in the bokeh column data source form
    SrcDict['set_1']=setMapDataSrc(Stats,'stlat','stlon',['type'])
    SrcDict['set_2']=setMapDataSrc(Eps,'eplat','eplon',['type'])

    # set plot attributes
    SrcDict['set_1_attributes']={'line_clr':'black','fill_clr':RGB(255,0,144),
                                'mrk':'inverted_triangle','size':8,'alpha':1}
    SrcDict['set_2_attributes']={'line_clr':'navy','fill_clr':'yellow',
                                'mrk':'circle_x','size':5,'alpha':0.6}

    return SrcDict

def setMapDataSrc(MapData,lat_field,lon_field,descriptions):
    datadict=dict()
    datadict['lon']=MapData[lon_field]
    datadict['lat']=MapData[lat_field]
    for desc in descriptions:
        datadict[desc]=MapData[desc]

    DataSrc = ColumnDataSource(data=datadict)

    return DataSrc

def setSWScatter3(pdDF):
    SrcDict=dict()
    # set plot attributes
    SrcDict['attrs']={'line_clr':RGB(0,0,200),'fill_clr':RGB(0,0,200),
                    'mrk':'circle','size':6,'alpha':0.5}
    SrcDict['datasrc']=setGenScatterSrc(pdDF,'delobs_1','delobsphase_1','delobs_2','delobsphase_2')
    return SrcDict

def setGenScatterSrc(pdDF,x_name,x_field,y_name,y_field,descriptions=None):
    datadict=dict()
    # print pdDF.head()
    # print x_name,x_field,y_name,y_field
    datadict[x_name]=pdDF[x_field]
    datadict[y_name]=pdDF[y_field]
    if descriptions is not None:
        for desc in descriptions:
            datadict[desc]=pdDF[desc]

    DataSrc = ColumnDataSource(data=datadict)
    return DataSrc

def getUserSelection_1(AppletType):
    bokehData=PreText(text='Select Summary or Raw Data (placeholder)', width=600)
    return bokehData

def getUserSelection_2(AppletType,DataSetList):
    if AppletType=='SurfaceWaves':
        bokehData=PreText(text="Selections:\n  Frequency\n  Study 1\n  Study 2\n\nPath:\n Epicenter Box\n Station Box", width=200)
    else:
        bokehData=PreText(text="Selection Widget", width=200)

    return bokehData

def getStatNetList(SW_stats):

    if SW_stats['stat_1'].str.contains('Null').any():
        # print("pulling stations from file 2")
        Nets=SW_stats['stat_2']
        Stats=SW_stats['stat_2']
    else:
        # print("pulling stations from file 1")
        Nets=SW_stats['stat_1']
        Stats=SW_stats['stat_1']

    # identify the format
    netstatsplits=['_','-']
    splitchar='_' # default
    for sp in netstatsplits:
        if Nets.str.contains(sp).any():
            splitchar=sp

    # have to do this in a loop because the data isn't well normalized
    NetList=[]
    StatList=[]
    for value in Nets:
        while '--' in value:
            value=value.replace('--','-')
        while '__' in value:
            value=value.replace('__','_')
        if splitchar in value:
            net=value.split(splitchar)[1]
            stat=value.split(splitchar)[0]
        else: # assume final 3 chars are the network
            net=value[-3:-1]
            stat=value[0:-4]
        NetList=NetList+[net]
        StatList=StatList+[net+'_'+stat]

    # sort the lists and add 'All' option.
    NetList=list(set(NetList))
    NetList.sort()
    NetList=['All']+NetList
    StatList.sort()
    StatList=['All']+StatList

    return StatList,NetList

def populateUserSelection(DataSetList,SW_stats):
    opts={}
    # studies
    opts['studies']=pd.unique(DataSetList['study']).tolist() # available studies

    # things that (COULD) depend on the studies chosen
    opts['StatList'],opts['NetList']=getStatNetList(SW_stats) # station list
    opts['path_R']=pd.unique(DataSetList['path']).tolist() # R1, R2, etc
    opts['normalmode']=[str(int(i)) for i in pd.unique(DataSetList['normalmode']).tolist()] # normal mode
    opts['period_min']=DataSetList['period_s'].min()
    opts['period_max']=DataSetList['period_s'].max()
    opts['period_step']=50

    return opts

def populateDefaultValues():
    ''' populate default studies for initial calculation. this could be smarter '''
    defs={}
    fi1='Lyon15.0.R1.100s.interp.rem3d'
    fi2='MBS11.0.R1.100s.interp.rem3d'

    defs['study1_file']=fi1
    defs['study2_file']=fi2
    defs['study1_group']=fi1.split('.')[0]
    defs['study2_group']=fi2.split('.')[0]
    period=fi1.split('.')[3]
    defs['period']=int(period.replace('s',''))
    defs['path_R']=fi1.split('.')[2]
    defs['normalmode']=str(int(fi1.split('.')[1]))

    return defs

def findFile(DataSetList,study,period,normalmode,path_R):

    matched=DataSetList.file_name[(DataSetList.study==study)& (DataSetList.period_s==period)
        & (DataSetList.path==path_R) &  (DataSetList.normalmode==int(normalmode))].tolist()

    if len(matched)>1:
        print("Warning, multiple matches found:")
        print(matched)
        print("selecting first")
        filename=matched[0]
    elif len(matched)<1:
        ficheck=study+"."+normalmode+'.'+path_R+'.'+str(period)+".interp.rem3d"
        print("Warning, NO file matches found for "+ ficheck)
        filename=''
    else:
        filename=matched[0]
        print("found file: "+filename)

    print 'filename is'
    print filename
    return filename

    # def updateSelected(self,inds):
    #     if len(inds) == 0:
    #         Nstats=len(self.Data)
    #         selectednets = list(set(self.Data['data_id'])) # pulls unique networks
    #     else:
    #         Nstats=len(inds)
    #         selectednets = list(set(self.Data['data_id'])) # pulls unique networks
    #
    #     selectednets.sort()
    #     stattext="Selected Station Details\n\nStation Count: "+str(Nstats)+"\nReferences Selected:"
    #     for net in selectednets:
    #         stattext=stattext+"\n    "+net
    #     stattext=stattext+"\n\nSlider: depth\nSlider:?\nButton: Calculate\nButton: Download"
    #     return stattext;

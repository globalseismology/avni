''' Surface Waves Applet '''
# import system classes
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer, HoverTool, ResetTool,BoxZoomTool
from bokeh.plotting import figure, curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Slider, Select, TextInput,MultiSelect,Button
import pandas as pd
import numpy as np

# import custom classes
from data.rem3d_wrapper import data_fetch
import data.processBokehData as pbd

# default data pull
Data=data_fetch('SurfaceWaves')

# set the files we want for default plot
defs=pbd.populateDefaultValues()
files_to_plot={1:defs['study1_file'],2:defs['study2_file']}

# pull those files into memory and process as needed
for fi in files_to_plot:
    Data.pullData(files_to_plot[fi],fi,False)
Data.processData() # prep the data for plotting

print 'rendering plots'
# creat map/station plot
TOOLS=[BoxZoomTool(),ResetTool()]
StationPlot = figure(tools=TOOLS, plot_width=600, plot_height=400, min_border=10, min_border_left=50,
           toolbar_location="above", x_axis_location=None, y_axis_location=None,
           title=Data.MapData['title'])
StationPlot.background_fill_color = "rgb(240,240,240)"

# plot first set (stations)
set='set_1'
attrs=set+'_attributes'
data_source1=Data.MapData[set]
rStationPlot1=StationPlot.scatter('lon', 'lat',source=data_source1,
        marker=Data.MapData[attrs]['mrk'], line_color=Data.MapData[attrs]['line_clr'],
        fill_color=Data.MapData[attrs]['fill_clr'], size=Data.MapData[attrs]['size'],
        alpha=Data.MapData[attrs]['alpha'])

 # plto second set (epicenters)
set='set_2'
attrs=set+'_attributes'
data_source2=Data.MapData[set]
rStationPlot2=StationPlot.scatter('lon', 'lat',source=data_source2,
        marker=Data.MapData[attrs]['mrk'], line_color=Data.MapData[attrs]['line_clr'],
        fill_color=Data.MapData[attrs]['fill_clr'], size=Data.MapData[attrs]['size'],
        alpha=Data.MapData[attrs]['alpha'])

# scatter plot 1
Scatter1=figure(tools=[],plot_width=300, plot_height=300, min_border=10, min_border_left=50,
                 title='Group 1')

# scatter plot 2
Scatter2=figure(tools=[],plot_width=300, plot_height=300, min_border=10, min_border_left=50,
                 title='Group 2')

# scatter plot 3
Scatter3=figure(tools=[],plot_width=300, plot_height=300, min_border=10, min_border_left=50,
                 title='Group 2 vs Group 1')
Scat3DataSource=Data.Scatter3Src['datasrc']
Attrs=Data.Scatter3Src['attrs']
rScatter3=Scatter3.scatter('delobs_1','delobs_2',source=Scat3DataSource,
        marker=Attrs['mrk'], line_color=Attrs['line_clr'],
        fill_color=Attrs['fill_clr'], size=Attrs['size'],
        alpha=Attrs['alpha'])

# create user selection panel (raw/summary)
rUserSelection_1 = pbd.getUserSelection_1(Data.AppletType)

# create user selection panel (studies, period, network, station)
opts=pbd.populateUserSelection(Data.DataSetList,Data.SW_stats)
period = Slider(title="Period [sec]", value=defs['period'], start=opts['period_min'],
                end=opts['period_max'], step=opts['period_step'])
study_group1 = Select(title="Group 1", value=defs['study1_group'],options=opts['studies'])
study_group2 = Select(title="Group 2", value=defs['study2_group'],options=opts['studies'])
mode = Select(title="Normal Mode", value=defs['normalmode'],options=opts['normalmode'])
path_R = Select(title="Path", value=defs['path_R'],options=opts['path_R'])
network= MultiSelect(title="Network (placeholder)", value=['All'],options=opts['NetList'],size=3)
stations = MultiSelect(title="Stations (placeholder)", value=['All'],options=opts['StatList'],size=3)

submit=Button(label='Submit',button_type='success')

sizing_mode = 'fixed'  # 'scale_width'
rUserSelection_2 = widgetbox([study_group1,study_group2,mode,path_R,period,submit,network,stations], sizing_mode=sizing_mode)

def select_values():
    # pull selected values
    network_vals= network.value
    station_vals= stations.value
    study1=str(study_group1.value)
    study2=str(study_group2.value)
    per=int(period.value)
    modeval=str(mode.value)
    pathval=path_R.value

    # find the corresponding file
    file_1=pbd.findFile(Data.DataSetList,study1,per,modeval,pathval)
    file_2=pbd.findFile(Data.DataSetList,study2,per,modeval,pathval)

    # put into dictionary
    selected={'period':per,'study1':study1,'study2':study2,
              'network':network_vals,'station':station_vals,'file_1':file_1,'file_2':file_2}
    return selected

def update():
    selected_vals=select_values() # get the selected values

    # pull the data files and process
    Data=data_fetch('SurfaceWaves')
    Data.pullData(selected_vals['file_1'],1,False)
    Data.pullData(selected_vals['file_2'],2,False)
    Data.processData() # will update DataSetList, Data.MapData, Data.ScatterData

    # update data sources
    data_source1.data=Data.MapData['set_1'].data
    data_source2.data=Data.MapData['set_2'].data
    Scat3DataSource.data=Data.Scatter3Src['datasrc'].data

    # update station & network list if files have changed
    StatList,NetList=pbd.getStatNetList(Data.SW_stats)
    network.options=NetList
    stations.options=StatList

    print("finished reprocessing")
    return

submit.on_click(update)

l=layout(
    [
    [rUserSelection_1],
    [StationPlot,rUserSelection_2],
    [Scatter1,Scatter2,Scatter3]
    ]
    )

curdoc().add_root(l)
curdoc().title = "REM3D"

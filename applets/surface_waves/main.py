''' Surface Waves Applet '''
# import system classes
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer, HoverTool, ResetTool
from bokeh.plotting import figure, curdoc, ColumnDataSource
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import PreText
import numpy as np

# import custom classes
from data.rem3d_wrapper import data_fetch

# initialize surface wave data object
Data=data_fetch('SurfaceWaves')
Data.rem3ddir='/home/chavlin/src/rem3d/rem3d' # issue with path definition!

# set the files we want to plot (this would come from user input)
files_to_plot={'file_1':'GDM52.0.R1.100s.interp.rem3d',
               'file_2':'MBS11.0.R1.100s.interp.rem3d'}

# pull those files into memory
for fi in files_to_plot:
    Data.pullData(files_to_plot[fi],fi,False)

# creat map/station plot
Data.buildMapData()
StatSrc = ColumnDataSource(data=dict(
    lon=Data.MapData['stlon'], lat=Data.MapData['stlat'], names=Data.MapData['stat'],
    datasource=Data.MapData['data_id'],
))
hover = HoverTool(tooltips=[ ("index", "$index"), ("(lat,lon)", "(@lat, @lon)"),
                     ("Station", "@names"), ("Data Source", "@datasource"),
])
resetbutton=ResetTool() # refresh doesn't call update (bug?)
TOOLS=[hover,BoxSelectTool(), LassoSelectTool(),resetbutton]

StationPlot = figure(tools=[], plot_width=800, plot_height=400, min_border=10, min_border_left=50,
           toolbar_location="above", x_axis_location=None, y_axis_location=None,
           title="Station Map")
StationPlot.background_fill_color = "rgb(240,240,240)"
StationPlot.select(BoxSelectTool).select_every_mousemove = False
StationPlot.select(LassoSelectTool).select_every_mousemove = False
rStationPlot=StationPlot.circle('lon','lat',size=8,source=StatSrc,selection_color='firebrick',selection_line_color='firebrick',
                           nonselection_fill_color='blue',nonselection_line_color="blue")

# 2. Station List div. Collect data from select stations, display.
rStationList = PreText(text='Selected Station Details', width=200)

# 3. file 1 renderer
Scatter1Src = ColumnDataSource(data=dict(
    distkm=Data.Data['distkm'][np.where(Data.Data['data_id']=='file_1')],
    delobsphase=Data.Data['delobsphase'][np.where(Data.Data['data_id']=='file_1')]
))
Scatter1=figure(tools=[],plot_width=400, plot_height=400, min_border=10, min_border_left=50,
                title='delobsphase vs distkm for '+files_to_plot['file_1'])
rScatter1=Scatter1.scatter('distkm','delobsphase',source=Scatter1Src)

# 4. file 2 renderer
Scatter2Src = ColumnDataSource(data=dict(
    distkm=Data.Data['distkm'][np.where(Data.Data['data_id']=='file_2')],
    delobsphase=Data.Data['delobsphase'][np.where(Data.Data['data_id']=='file_2')]
))
Scatter2=figure(tools=[],plot_width=400, plot_height=400, min_border=10, min_border_left=50,
                title='delobsphase vs distkm for '+files_to_plot['file_2'])
rScatter2=Scatter2.scatter('distkm','delobsphase',source=Scatter2Src)


# Layout of all the plots
l=layout([[StationPlot],[Scatter1,Scatter2]])
curdoc().add_root(l)
curdoc().title = "REM3D"

# handle requests
# def update_selected_stations(attr, old, new):
#     inds=np.array(new['1d']['indices'])
#     rStationList.text=Data.updateSelected(inds)
#
# rStationPlot.data_source.on_change('selected', update_selected_stations)

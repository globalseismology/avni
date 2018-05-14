'''
minimal working example of a bokeh webapp using subdirectories only
'''
# import standard classes
import numpy as np
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer, HoverTool, ResetTool
from bokeh.plotting import figure, curdoc, ColumnDataSource
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import PreText

# import custom classes
from data.data_fetch import data_fetch

# Data Prep: pull station data
stations=data_fetch()
stations.set_station_latlons() # lat,lon,station_name,station_network

# 1. StationPlot. Rollover & select scatter plot of station
StatSrc = ColumnDataSource(data=dict(
    lon=stations.lon, lat=stations.lat, names=stations.station_name,
    networks=stations.station_network, netcolor=stations.color_map,
))
hover = HoverTool(tooltips=[ ("index", "$index"), ("(lat,lon)", "(@lat, @lon)"),
                     ("Station ID", "@names"), ("Station Network", "@networks"),
])
resetbutton=ResetTool() # refresh doesn't call update (bug?)
TOOLS=[hover,BoxSelectTool(), LassoSelectTool(),resetbutton]

StationPlot = figure(tools=TOOLS, plot_width=600, plot_height=600, min_border=10, min_border_left=50,
           toolbar_location="above", x_axis_location=None, y_axis_location=None,
           title="Station Map")
StationPlot.background_fill_color = "rgb(240,240,240)"
StationPlot.select(BoxSelectTool).select_every_mousemove = False
StationPlot.select(LassoSelectTool).select_every_mousemove = False
rStationPlot=StationPlot.circle('lon','lat',size=10,source=StatSrc,color='netcolor')

# 2. Station List div. Collect data from select stations, display.
rStationList = PreText(text='Selected Station Details', width=200)

# Layout of all the plots
l=layout([[rStationList,StationPlot]])
curdoc().add_root(l)
curdoc().title = "REM3D Explorer"

# handle requests
def update_selected_stations(attr, old, new):
    inds=np.array(new['1d']['indices'])
    rStationList.text=stations.update_selected(inds)

rStationPlot.data_source.on_change('selected', update_selected_stations)

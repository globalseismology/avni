import os
import sys
import imp
import cPickle
import numpy as np
import pyqtgraph as pg
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
from PyQt4 import uic
from PyQt4 import QtGui, QtCore
from PyQt4.QtGui import QVBoxLayout, QFileDialog
from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
from matplotlib import gridspec

from rem3d.geolib import ddelazgc
from rem3d.mapping import spher2cart
from rem3d.plots import globalmap
from rem3d.plots import getmodeltransect
from rem3d.plots import setup_axes

DIR = os.path.dirname(os.path.abspath(__file__))
dbs_path = '/home/romaguir/dbs'

class Window(QtGui.QMainWindow):
    
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = layout.Ui_MainWindow()
        self.ui.setupUi(self)
        self.k = 8
        #self.nc4_model = xr.open_dataarray('~/dbs/S40RTS_test.nc4')
        self.nc4_model = xr.open_dataarray('~/dbs/S40RTS_pixel_0.5x0.5.nc4')
        print 'READING TREE'
        self.kdtree = cPickle.load(open(dbs_path+'/S40RTS_pixel_0.5x0.5.KDTree.pkl'))
        print 'DONE READING TREE'
        self.global_map_figure = self.ui.global_map.fig
        self.global_map_ax = self.global_map_figure.add_axes([0.01,0.01,0.99,0.99])
        self.bmap_globe = None
        self.region_map_figure = self.ui.local_map.fig
        self.region_map_ax = self.region_map_figure.add_axes([0.01,0.01,0.99,0.99])
        self.bmap_region = None
        self.cross_section_figure = self.ui.cross_section.fig
        self.cross_section_ax = self.cross_section_figure.add_axes([0.01,0.01,0.99,0.99])

        #cross section endpoints
        self.xsec_lat1 = None
        self.xsec_lat2 = None
        self.xsec_lon1 = None
        self.xsec_lon2 = None
        self.endpoint1 = None
        self.endpoint2 = None
        self.gcpath_r = None

        #initialize plots
        self.plot_global_map()

        #setup click events
        cid = self.ui.local_map.mpl_connect('button_press_event', self.onclick)

    def plot_global_map(self):
        plt.cla()
        dv_array = self.nc4_model.interp(depth=float(self.ui.depth_slider.value()))
        #dv_array = self.nc4_model.interp(depth=100.0)
        if self.bmap_globe == None:
            self.bmap_globe = Basemap(ax = self.global_map_ax,projection='robin',
                        lat_0=0,lon_0=0,resolution='c')
            self.bmap_globe.drawcoastlines(zorder=90)

        self.bmap_globe.imshow(dv_array.T,vmin=-2.0,vmax=2.0,cmap='jet_r')
        self.global_map_figure.canvas.draw()
        self.draw_screen_poly()
        self.global_map_figure.canvas.flush_events()

    def plot_regional_map(self):
        plt.cla()
        if self.bmap_region == None:
            self.bmap_region = Basemap(ax = self.region_map_ax,projection='merc',
                                       llcrnrlat=self.ui.lat_min_slider.value(),
                                       llcrnrlon=self.ui.lon_min_slider.value(),
                                       urcrnrlat=self.ui.lat_max_slider.value(),
                                       urcrnrlon=self.ui.lon_max_slider.value(),
                                       resolution='c') 
            meridians = np.arange(-180,181,10)
            self.bmap_region.drawmeridians(meridians,labels=[True,False,False,False])
            parallels = np.arange(-90,91,10)
            self.bmap_region.drawparallels(parallels,labels=[True,False,False,False])

        self.bmap_region.drawcoastlines(zorder=90)
        self.region_map_figure.canvas.draw()
        self.region_map_figure.canvas.flush_events()

    def draw_screen_poly(self):
        plt.cla()
        lon_min = self.ui.lon_min_slider.value()
        lon_max = self.ui.lon_max_slider.value()
        lat_min = self.ui.lat_min_slider.value()
        lat_max = self.ui.lat_max_slider.value()

        x1, y1 = self.bmap_globe([lon_min,lon_max],[lat_min,lat_min] )
        x2, y2 = self.bmap_globe([lon_min,lon_max],[lat_max,lat_max] )
        line1 = self.bmap_globe.drawgreatcircle(lon_min,lat_min,lon_min,lat_max,color='purple',linewidth=3,zorder=99)
        line2 = self.bmap_globe.drawgreatcircle(lon_max,lat_min,lon_max,lat_max,color='purple',linewidth=3,zorder=99)
        line3 = self.bmap_globe.plot(x1,y1,color='purple',linewidth=3,zorder=99)
        line4 = self.bmap_globe.plot(x2,y2,color='purple',linewidth=3,zorder=99)

        self.global_map_figure.canvas.draw()
        line1[0].remove()
        line2[0].remove()
        line3[0].remove()
        line4[0].remove()
        self.global_map_figure.canvas.flush_events()

    def on_depth_slider_valueChanged(self):
        print self.ui.depth_slider.value()
        self.plot_global_map()

    def on_lon_min_slider_valueChanged(self):
        self.draw_screen_poly()

    def on_lat_min_slider_valueChanged(self):
        self.draw_screen_poly()

    def on_lon_max_slider_valueChanged(self):
        self.draw_screen_poly()

    def on_lat_max_slider_valueChanged(self):
        self.draw_screen_poly()

    def on_submit_button_released(self):
 
        self.region_map_ax.clear()
        self.region_map_ax = self.region_map_figure.add_axes([0.01,0.01,0.99,0.99])
        self.bmap_region = None
        self.plot_regional_map()

    def onclick(self,event):
        #print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
        #      ('double' if event.dblclick else 'single', event.button,
        #       event.x, event.y, event.xdata, event.ydata))
        xpt = event.xdata
        ypt = event.ydata
        lonpt, latpt = self.bmap_region(xpt,ypt,inverse=True)

        if self.xsec_lat1 == None:
            self.xsec_lat1 = latpt
            self.xsec_lon1 = lonpt
            self.endpoint1 = self.bmap_region.scatter(xpt,ypt,c='b',marker='o')

        elif self.xsec_lat1 != None and self.xsec_lat2 == None:
            self.xsec_lat2 = latpt
            self.xsec_lon2 = lonpt
            self.endpoint2 = self.bmap_region.scatter(xpt,ypt,c='b',marker='o')
            self.gcpath_r = self.bmap_region.drawgreatcircle(self.xsec_lon1,
                                                             self.xsec_lat1,
                                                             self.xsec_lon2,
                                                             self.xsec_lat2,
                                                             c='b',linewidth=2,zorder=99)

        elif self.xsec_lat1 != None and self.xsec_lat2 != None:
            self.xsec_lat1 = None
            self.xsec_lon1 = None
            self.xsec_lat2 = None
            self.xsec_lon2 = None
            self.endpoint1.remove()
            self.endpoint2.remove()
            self.gcpath_r[0].remove()

        print 'POINT 1 = ', self.xsec_lon1, self.xsec_lat1
        print 'POINT 2 = ', self.xsec_lon2, self.xsec_lat2

        self.region_map_figure.canvas.draw()
        self.region_map_figure.canvas.flush_events()

    def on_cross_section_button_released(self):

        self.cross_section_ax.clear()
        gs = gridspec.GridSpec(1,1)
        numevalx = 200
        numevalz = 200
        #self.cross_section_ax = self.cross_section_figure.add_axes([0.01,0.01,0.99,0.99])
        gcdelta,az1,az2 = ddelazgc(self.xsec_lat1,self.xsec_lon1,
                                   self.xsec_lat2,self.xsec_lon2)
        evalpoints,xsec,junk = getmodeltransect(lat1=self.xsec_lat1,lng1=self.xsec_lon1,
                               azimuth=az1,gcdelta=gcdelta,f=self.nc4_model,
                               tree=self.kdtree,numevalx=numevalx,numevalz=numevalz)

        #theta = [0.0,gcdelta]
        theta = [90.0-gcdelta/2., 90.+gcdelta/2.]
        radius = [3480.0, 6371.0]
        self.cross_section_ax, aux_ax = setup_axes(self.cross_section_figure, gs[0], theta,
                                                   radius = radius ,numdegticks=10)
        self.cross_section_ax.set_aspect('equal')
        aux_ax.set_aspect('equal')
        #interp_values = np.array(xsec) #maybe transpose?
        #interp_values = (interp_values - interp_values.mean()).reshape(numevalx,numevalz)
        #self.cross_section_ax.imshow(np.flipud(xsec),vmin=-2.0,vmax=2.0,aspect='auto',cmap='jet_r')
        grid_x, grid_y = np.meshgrid(np.linspace(theta[0],theta[1],numevalx),np.linspace(radius[0],radius[1],numevalz))
        #self.cross_section_ax.pcolormesh(grid_x,grid_y,interp_values,vmin=-2.0,vmax=2.0,cmap='jet_r')
        aux_ax.pcolormesh(grid_x,grid_y,xsec,vmin=-2.0,vmax=2.0,cmap='jet_r')
        self.cross_section_figure.canvas.draw()
        self.cross_section_figure.canvas.flush_events()

def use_ui_layout():
    DIR = os.path.dirname(os.path.abspath(__file__)) 
    print DIR
    ui_file = DIR+'/layout.ui'
    py_file = DIR+'/layout.py'
    with open(py_file, 'w') as open_file:
        uic.compileUi(ui_file, open_file)

    import_name = os.path.splitext(os.path.basename(py_file))[0]
    print "import name", import_name
    globals()[import_name] = imp.load_source(import_name, py_file)

def launch():
    
    use_ui_layout()
    app = QtGui.QApplication(sys.argv, QtGui.QApplication.GuiClient)
    window = Window()
    window.show()
    app.installEventFilter(window)
    window.raise_()                                                                                                                                                                   
    os._exit(app.exec_())

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

from rem3d.f2py import ddelazgc
from rem3d.mapping import spher2cart
from rem3d.plots import globalmap
from rem3d.plots import standardcolorpalette
from rem3d.plots import getmodeltransect
from rem3d.plots import setup_axes
from rem3d.plots import plot_hotspots
from rem3d.plots import plot_plates
from rem3d import tools

DIR = os.path.dirname(os.path.abspath(__file__))
dbs_path = '/home/romaguir/dbs'
plt.box('False')

class Window(QtGui.QMainWindow):
    
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = layout.Ui_MainWindow()
        self.ui.setupUi(self)
        self.k = 8
        #self.nc4_model = xr.open_dataarray('~/dbs/S40RTS_test.nc4')
        #self.nc4_model = xr.open_dataarray('~/dbs/S40RTS_pixel_0.5x0.5.nc4')
        self.nc4_model = xr.open_dataarray('~/dbs/S362ANI+M_vs.nc4')
        print 'READING TREE'
        self.kdtree = cPickle.load(open(dbs_path+'/S362ANI+M.KDTree.pkl'))
        print 'DONE READING TREE'
        self.global_map_figure = self.ui.global_map.fig
        self.global_map_ax = self.global_map_figure.add_axes([0.05,0.1,0.75,0.8])
        self.bmap_globe = None
        self.region_map_figure = self.ui.local_map.fig
        self.region_map_ax = self.region_map_figure.add_axes([0.1,0.1,0.8,0.8])
        self.region_map_ax.set_visible(False)
        self.bmap_region = None
        self.cross_section_figure = self.ui.cross_section.fig
        self.cross_section_ax = self.cross_section_figure.add_axes([0.1,0.1,0.8,0.8])
        self.cross_section_ax.set_visible(False)
        self.global_map_cbar = None
        self.cbar_ax = None
        self.cross_section_drawn = False

        self.hotspots_plotted = False
        self.coastlines_plotted = False
        self.plates_plotted = False

        self.hotspots = None
        self.coastlines = None
        self.plates = None
        self.ridge = None
        self.transform = None
        self.trench = None

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
            #self.bmap_globe.drawcoastlines(color='darkgray')
            if self.coastlines_plotted == False:
                self.coastlines_plotted = True
                print 'SELF.COASTLINES SHOULD BE DEFINED'
                self.coastlines = self.bmap_globe.drawcoastlines(color=(0.25,0.25,0.25))

        rem3d_cpt = standardcolorpalette(name='rem3d')
        print rem3d_cpt

        #self.bmap_globe = globalmap(self.global_map_ax,dv_array,vmin=-2.0,vmax=2.0,dbs_path=dbs_path)
        vmax = self.ui.cbar_max.value()
        vmin = -vmax
   
        im = self.bmap_globe.imshow(np.roll(dv_array.T,180),vmin=vmin,vmax=vmax,cmap=rem3d_cpt)
        #im = self.bmap_globe.contourf(np.roll(dv_array.T,180),vmin=vmin,vmax=vmax,cmap=rem3d_cpt)
        self.global_map_ax.set_title('depth = {} km'.format(self.ui.depth_slider.value()))

        if self.global_map_cbar == None:
        #    #self.global_map_cbar = self.global_map_figure.colorbar(im, ax=self.global_map_ax,label='$\delta V_S / V_S$ %')
        #    self.cbar_ax = self.global_map_figure.add_axes([0.85,0.1,0.05,0.8])
            self.cbar_ax = self.global_map_figure.add_axes([0.825,0.1,0.025,0.8])
            self.global_map_cbar = self.global_map_figure.colorbar(im, cax=self.cbar_ax,label='$\delta V_S / V_S$ %')
        else:
            self.global_map_cbar.remove() 
            self.cbar_ax = self.global_map_figure.add_axes([0.825,0.1,0.025,0.8])
            self.global_map_cbar = self.global_map_figure.colorbar(im, cax=self.cbar_ax,label='$\delta V_S / V_S$ %')

        self.global_map_figure.canvas.draw()
        self.draw_screen_poly()
        self.global_map_figure.canvas.flush_events()
        im = None
        #self.global_map_cbar.remove()

    def plot_regional_map(self):
        plt.cla()
        dv_array = self.nc4_model.interp(depth=float(self.ui.depth_slider.value()))
        if self.bmap_region == None:
            self.bmap_region = Basemap(ax = self.region_map_ax,projection='merc',
                                       llcrnrlat=self.ui.lat_min_slider.value(),
                                       llcrnrlon=self.ui.lon_min_slider.value(),
                                       urcrnrlat=self.ui.lat_max_slider.value(),
                                       urcrnrlon=self.ui.lon_max_slider.value(),
                                       resolution='c') 

            lon_width = self.ui.lon_max_slider.value() - self.ui.lon_min_slider.value()

            if lon_width >= 180:
                meridians = np.arange(-180,181,30)
                parallels = np.arange(-90,91,30)
            elif lon_width < 180.0 and lon_width > 90.0:
                meridians = np.arange(-180,181,20)
                parallels = np.arange(-90,91,20)
            else:
                meridians = np.arange(-180,181,10)
                parallels = np.arange(-90,91,10)

            self.bmap_region.drawmeridians(meridians,labels=[True,False,True,False])
            self.bmap_region.drawparallels(parallels,labels=[True,False,True,False])

        #to show velocity model
        #rem3d_cpt = standardcolorpalette(name='rem3d')
        #self.bmap_region.imshow(dv_array,cmap=rem3d_cpt,vmin=-2.0,vmax=2.0)

        self.bmap_region.drawcoastlines(color='darkgray')
        self.region_map_ax.set_visible(True)
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
        line1 = self.bmap_globe.drawgreatcircle(lon_min,lat_min,lon_min,lat_max,color='magenta',linewidth=3,zorder=99)
        line2 = self.bmap_globe.drawgreatcircle(lon_max,lat_min,lon_max,lat_max,color='magenta',linewidth=3,zorder=99)
        line3 = self.bmap_globe.plot(x1,y1,color='magenta',linewidth=3,zorder=99)
        line4 = self.bmap_globe.plot(x2,y2,color='magenta',linewidth=3,zorder=99)

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
        self.region_map_ax = self.region_map_figure.add_axes([0.1,0.1,0.8,0.8])
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
            self.endpoint1 = self.bmap_region.scatter(xpt,ypt,marker='o',edgecolor='black',facecolor='white',zorder=90)

        elif self.xsec_lat1 != None and self.xsec_lat2 == None:
            self.xsec_lat2 = latpt
            self.xsec_lon2 = lonpt
            self.endpoint2 = self.bmap_region.scatter(xpt,ypt,color='magenta',marker='o',zorder=90)
            self.gcpath_r = self.bmap_region.drawgreatcircle(self.xsec_lon1,
                                                             self.xsec_lat1,
                                                             self.xsec_lon2,
                                                             self.xsec_lat2,
                                                             color='cyan',linewidth=2,zorder=89)

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

    #def on_cross_section_button_released(self):
    def draw_cross_section(self):

        self.cross_section_ax.clear()
        gs = gridspec.GridSpec(1,1)
        numevalx = 200
        numevalz = 200
        #self.cross_section_ax = self.cross_section_figure.add_axes([0.01,0.01,0.99,0.99])
        gcdelta,az1,az2 = ddelazgc(self.xsec_lat1,self.xsec_lon1,
                                   self.xsec_lat2,self.xsec_lon2)
        evalpoints,xsec,junk = getmodeltransect(lat1=self.xsec_lat1,lng1=self.xsec_lon1,
                               azimuth=az1,gcdelta=gcdelta,f=self.nc4_model,
                               tree=self.kdtree,numevalx=numevalx,numevalz=numevalz,
                               k=self.k)

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
        rem3d_cpt = standardcolorpalette(name='rem3d')
        vmax = self.ui.cbar_max.value()
        vmin = -vmax
        aux_ax.pcolormesh(grid_x,grid_y,np.fliplr(xsec),vmin=vmin,vmax=vmax,cmap=rem3d_cpt)
        #self.cross_section_figure.tight_layout()
        self.cross_section_drawn = True
        self.cross_section_figure.canvas.draw()
        self.cross_section_figure.canvas.flush_events()


    def on_cross_section_button_released(self):
        self.draw_cross_section()

    def on_cbar_max_valueChanged(self):
        #self.global_map_cbar = None
        dv_array = self.nc4_model.interp(depth=float(self.ui.depth_slider.value()))
        rem3d_cpt = standardcolorpalette(name='rem3d')
        vmin=-self.ui.cbar_max.value()
        vmax=self.ui.cbar_max.value()
        im = self.bmap_globe.imshow(np.roll(dv_array.T,180),vmin=vmin,vmax=vmax,cmap=rem3d_cpt)
        #self.global_map_figure.canvas.draw()
        #self.global_map_figure.canvas.flush_events()
        self.plot_global_map()
        #self.draw_screen_poly()

        #update cross section (if it exists)
        if self.cross_section_drawn:
            self.draw_cross_section()

    def on_hotspots_box_stateChanged(self):
        if self.hotspots_plotted == False:
            self.hotspots_plotted = True
            #plot_hotspots(self.bmap_globe, dbs_path = dbs_path)
            hotspots = tools.readjson('%s/hotspots.json' % (dbs_path))
            x, y = self.bmap_globe(hotspots[:,0], hotspots[:,1])
            self.hotspots = self.bmap_globe.scatter(x, y,marker='^',facecolor='white',edgecolor='black')
        else:
            self.hotspots_plotted = False
            self.hotspots.remove()

        self.global_map_figure.canvas.draw()
        self.draw_screen_poly()
        self.global_map_figure.canvas.flush_events()

    def on_plate_boundaries_box_stateChanged(self):
        plot_plates(self.bmap_globe, dbs_path = dbs_path)

        #if self.plates_plotted == False:
        #    self.plate_plotted = True
        #    boundtypes=['ridge', 'transform', 'trench']
        #    for bound in boundtypes:
        #        name, segs = tools.readjson('%s/%s.json' % (dbs_path,bound))
        #        segs = np.array(segs)
        #        ind_nan = np.nonzero(np.isnan(segs[:,0]))
        #        segs[ind_nan,0] = 0
        #        segs[ind_nan,1] = 0
        #        x,y = self.bmap_globe(segs[:,0], segs[:,1])
        #        x[ind_nan] = np.nan
        #        y[ind_nan] = np.nan
        #        dx = np.abs(x[:1] - x[:-1])
        #        ind_jump, = np.nonzero(dx > 1.75e7)
        #        x[ind_jump] = np.nan

        #        if bound == 'ridge':
        #            self.ridge = self.bmap_globe.plot(x,y,color='red')
        #        elif bound == 'transform':
        #            self.transform = self.bmap_globe.plot(x,y,color='green')
        #        elif bound == 'trench':
        #            self.transform = self.bmap_globe.plot(x,y,color='blue')

            #xs = np.array(xs).flatten()
            #ys = np.array(ys).flatten()
            #self.plates = self.bmap_globe.plot(xs,ys,'-',color='cyan')

        #else:
        #    self.plates_plotted = False
            #self.plates.remove()

        self.global_map_figure.canvas.draw()
        self.draw_screen_poly()
        self.global_map_figure.canvas.flush_events()

    def on_coastlines_box_stateChanged(self):
        if self.coastlines_plotted == True:
            self.coastlines.remove()
            self.coastlines_plotted = False
        else:
            self.coastlines = self.bmap_globe.drawcoastlines(color=(0.25,0.25,0.25))
            self.coastlines_plotted = True

        self.global_map_figure.canvas.draw()
        self.draw_screen_poly()
        self.global_map_figure.canvas.flush_events()


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

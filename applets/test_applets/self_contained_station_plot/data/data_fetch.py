""" class for simple test case. fabricates some data for plotting """
import numpy as np

class data_fetch:

    def selection_histogram(self):
        # from the bokeh server example
        x1 = np.random.normal(loc=5.0, size=400) * 100
        y1 = np.random.normal(loc=10.0, size=400) * 10

        x2 = np.random.normal(loc=5.0, size=800) * 50
        y2 = np.random.normal(loc=5.0, size=800) * 10

        x3 = np.random.normal(loc=55.0, size=200) * 10
        y3 = np.random.normal(loc=4.0, size=200) * 10

        self.x = np.concatenate((x1, x2, x3))
        self.y = np.concatenate((y1, y2, y3))
        return;

    def set_station_latlons(self,station_networks='all'):
        # need to pull this data from somewhere, for now, let's invent some
        minlat=30
        maxlat=50
        minlon=260
        maxlon=270
        nstations=95
        self.lat=np.random.random_sample(np.zeros(nstations).shape)*(maxlat-minlat)+minlat
        self.lon=np.random.random_sample(self.lat.shape)*(maxlon-minlon)+minlon
        self.station_name=list()
        self.station_network=list()
        self.unique_networks=list()
        self.color_map=list() # a 3-tuple of integers (r,g,b) between 0 and 255
        color_sequence=["rgb(255,0,0)","rgb(255,0,255)","rgb(0,255,0)","rgb(0,0,255)","rgb(0,255,255)","rgb(0,0,0)"]
        current_color=color_sequence[0];
        statnum=1
        netnum=1
        colornum=0
        self.unique_networks.append('Network '+str(netnum))
        for stat in self.lat:
            self.station_name.append('Station '+str(statnum))
            self.station_network.append('Network '+str(netnum))
            self.color_map.append(current_color)
            statnum=statnum+1
            if statnum%20==0:
                netnum=netnum+1
                colornum=colornum+1
                if (colornum>len(color_sequence)-1):
                    colornum=0
                current_color=color_sequence[colornum]
                self.unique_networks.append('Network '+str(netnum))

        self.station_name=np.array(self.station_name)
        self.station_network=np.array(self.station_network)
        self.unique_networks=np.array(self.unique_networks)
        self.nnets=netnum
        self.nstats=len(self.lat)
        return;

    def update_selected(self,inds):
        if len(inds) == 0:
            Nstats=len(self.lat)
            selectednets = list(set(self.station_network)) # pulls unique networks
        else:
            Nstats=len(inds)
            selectednets = list(set(self.station_network[inds])) # pulls unique networks

        selectednets.sort()
        stattext="Selected Station Details\n\nStation Count: "+str(Nstats)+"\nNetworks Selected:"
        for net in selectednets:
            stattext=stattext+"\n    "+net
        stattext=stattext+"\n\nSlider: depth\nSlider:?\nButton: Calculate\nButton: Download"
        return stattext;

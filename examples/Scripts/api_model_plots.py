#!/usr/bin/env python
'''
example for using the avni api client: fetch model evaluations and plot a cross section

python api_model_explorer api_key

api_key is optional if api_initialization has been run.


'''

from avni.api.client import Client as r3d
from avni.api.model import Model
import matplotlib.pyplot as plt
import numpy as np
import argparse
from avni.plots.models import plot1section, plot1globalmap

def main():
    # parse
    parser = argparse.ArgumentParser(description='plots some profiles')
    parser.add_argument('-k', '--key', type=str,default='',help='api key')
    arg = parser.parse_args()

    # initialize API connection
    if arg.key=='':
        conn=r3d()
    else:
        conn=r3d(api_key=arg.key) # the connection object

    ModelInstance=Model(conn) # give the model instance the connection object

    # get the cross section data
    args={'lat':42.2,'lon':232.0,'azimuth':80.0,'gcdelta':85.}
    max_depth = 2900. 
    min_depth = 50. 
    args['radius_min_km']=6371-max_depth
    args['radius_max_km']=6371-min_depth
    args['quickInterp']=1
    xs=ModelInstance.crossSection(args)
                
    # now pass the cross section data to plot1section 
    xsec_data = xs['vs']
    vmin = -3 
    vmax = 3 
    plot1section(args['lat'],args['lon'],args['azimuth'],args['gcdelta'],None,'vs',vmin,vmax,xsec_data = xsec_data)

    # get the fixed depth values 
    fixed_d = ModelInstance.fixedDepth({'depth':250.})    

    # create the required named array 
    lons =fixed_d['lon']- 180. 
    long,latg=np.meshgrid(lons,fixed_d['lat'])
    data = np.vstack((latg.ravel(),long.ravel(),fixed_d['vs'].ravel())).transpose()
    dt = {'names':['latitude', 'longitude', 'value'], 'formats':[np.float, np.float, np.float]}
    valarray = np.zeros(len(data), dtype=dt)
    valarray['latitude'] = data[:,0]; valarray['longitude'] = data[:,1]; valarray['value'] = data[:,2]
    vmin = -5 
    vmax = 5
    plot1globalmap(valarray,vmin,vmax,colorpalette='avni')
if __name__== "__main__":
    main()

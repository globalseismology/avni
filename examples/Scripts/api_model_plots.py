#!/usr/bin/env python
'''
example for using the avni api client: fetch model evaluations and cross-sections

python api_model_explorer api_key

api_key is optional if api_initialization has been run.

(this example will include simple plotting for the call results)

'''

from avni.api.client import Client as r3d
from avni.api.model import Model
import matplotlib.pyplot as plt
import numpy as np
import argparse
from avni.plots.models import section

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
    args={'lat':42.2,'lon':232.0,'azimuth':80.0,'gcdelta':200.}
    max_depth = 2900. 
    min_depth = 50. 
    args['radius_min_km']=6371-max_depth
    args['radius_max_km']=6371-min_depth
    args['quickInterp']=1
    xs=ModelInstance.crossSection(args)
                
    # plot it by calling section directly
    fig = plt.figure()
    xsec_data = xs['vs']
    sec = section(fig,args['lat'],args['lon'],args['azimuth'],args['gcdelta'],None,'vs', xsec_data = xsec_data)
    plt.show() 
        
    # another 
    args={'lat':-20.2,'lon':32.0,'azimuth':32.0,'gcdelta':39.0}
    max_depth = 2900. 
    min_depth = 50. 
    args['radius_min_km']=6371-max_depth
    args['radius_max_km']=6371-min_depth
    args['quickInterp']=1
    xs=ModelInstance.crossSection(args)
    fig = plt.figure()
    xsec_data = xs['vs']
    sec = section(fig,args['lat'],args['lon'],args['azimuth'],args['gcdelta'],None,'vs', xsec_data = xsec_data)
    plt.show() 

if __name__== "__main__":
    main()

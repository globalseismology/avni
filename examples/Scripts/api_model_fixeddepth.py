#!/usr/bin/env python
'''
example for using the rem3d api client: fetch model cross-sections at fixed depths

python api_model_fixeddepth -k api_key

api_key is optional if api_initialization has been run.

(this example includes simple plotting for the call results)

'''

from rem3d.api.client import Client as r3d
from rem3d.api.model import Model
import numpy as np
import argparse
import matplotlib.pyplot as plt

def main():
    # parse
    parser = argparse.ArgumentParser(description='some model api calls')
    parser.add_argument('-k', '--key', type=str,default='',help='api key')
    arg = parser.parse_args()

    # initialize API connection
    if arg.key=='':
        conn=r3d()
    else:
        conn=r3d(api_key=arg.key) # the connection object

    ModelInstance=Model(conn) # give the model instance the connection object

    # evaluate model at a fixed depth, whole earth
    print("Fetching whole earth fixed depth")
    result=ModelInstance.fixedDepth({'depth':80.0})
    plt.figure()
    plt.subplot(1,3,1)
    plt.contourf(result['lon'],result['lat'],result['vs'],30)
    cR=[]
    cR.append(np.min(result['vs']))
    cR.append(np.max(result['vs']))
    plt.clim(cR)
    plt.title('80 km '+result['input_args']['model'])

    # limit the lat/lon ranges at fixed depth for N. America, adjust resolution
    args={
        'depth':80.0,'lat1':0.,'lat2':80,'lon1':-130,'lon2':-60,
        'Nlat':30,'Nlon':80
        }
    print("Fetching N. America fixed depth")
    result=ModelInstance.fixedDepth(args)
    plt.subplot(1,3,2)
    plt.contourf(result['lon'],result['lat'],result['vs'],30)
    plt.clim(cR)
    plt.title('North America, 80 km '+result['input_args']['model'])

    # use a different model
    # limit the lat/lon ranges at fixed depth for N. America, adjust resolution
    args['model']='S40RTS'
    print("Fetching N. America fixed depth, model 2")
    result2=ModelInstance.fixedDepth(args)
    plt.subplot(1,3,3)
    plt.contourf(result2['lon'],result2['lat'],result2['vs'],30)
    plt.clim(cR)
    plt.title('North America, 80 km '+result2['input_args']['model'])


    plt.show()


if __name__== "__main__":
    main()

#!/usr/bin/env python
'''
example for using the rem3d api client: fetch model evaluations and cross-sections

python api_model_explorer api_key

api_key is optional if api_initialization has been run.

(this example will include simple plotting for the call results)

'''

from rem3d.api.client import Client as r3d
from rem3d.api.model import Model
import numpy as np
import argparse

def main():
    # parse
    parser = argparse.ArgumentParser(description='f2py calls via API')
    parser.add_argument('-k', '--key', type=str,default='',help='api key')
    arg = parser.parse_args()

    # initialize API connection
    if arg.key=='':
        conn=r3d()
    else:
        conn=r3d(api_key=arg.key) # the connection object

    ModelInstance=Model(conn) # give the model instance the connection object

    # list available models & kernels
    model_list=ModelInstance.listModels()
    print("\nModel List")
    print(model_list)

    # evaluate model at a single point
    model_evaluation=ModelInstance.evaluate_points({'lat':42.2,'lon':232.0,'depth':80.0})
    print("\nModel Evaluation of scalar lat/lon/depth point:")
    print('vs='+str(model_evaluation['vs']))

    # evalualte at a list/array of points
    lat=np.linspace(45.0,60.0,10)
    dep=np.ones(lat.shape)*150
    lon=np.ones(lat.shape)*15
    args={'lat':lat,'lon':lon,'depth':dep}
    print("\nModel Evaluation of 1d numpy arrays")
    model_evaluation=ModelInstance.evaluate_points(args)
    print('vs='+str(model_evaluation['vs']))

    # get a cross-section

    # required args:
    gcdelta=100.0;
    args={'lat':42.2,'lon':232.0,'azimuth':80.0,'gcdelta':gcdelta}

    # optional args:
    args['radius_min_km']=6371-50.
    args['radius_max_km']=6371-200.

    print("\nquick interp crossection")
    args['quickInterp']=1
    xs=ModelInstance.crossSection(args)

    print("\nfull model load crossection")
    args['quickInterp']=0
    xs2=ModelInstance.crossSection(args)

    print("\nquick interp crossection with model specified")
    args['quickInterp']=1
    args['model']='S40RTS'
    xs=ModelInstance.crossSection(args)


if __name__== "__main__":
    main()

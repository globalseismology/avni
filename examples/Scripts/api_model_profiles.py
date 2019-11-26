#!/usr/bin/env python
'''
example for using the rem3d api client: fetch model evaluations and cross-sections

python api_model_explorer api_key

api_key is optional if api_initialization has been run.

(this example will include simple plotting for the call results)

'''

from rem3d.api.client import Client as r3d
from rem3d.api.model import Model
import matplotlib.pyplot as plt
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

    # get some depth profiles
    print("\ndepth profile with defaults")
    args={'lat':42.2,'lon':232.0}
    depth_profile=ModelInstance.depthProfile(args)

    print("\ndepth profile specified depth range")
    args={'lat':42.2,'lon':232.0,'depthMin':50,'depthMax':250,'N_depth':10}
    depth_profile_2=ModelInstance.depthProfile(args)

    ax=plt.axes()
    ax.plot(depth_profile['vs'],depth_profile['depth'],'k')
    ax.plot(depth_profile_2['vs'],depth_profile_2['depth'],'.r')
    ax.invert_yaxis()
    ax.set_xlabel('V_s perturbation')
    ax.set_ylabel('depth')
    plt.show()


if __name__== "__main__":
    main()

#!/usr/bin/env python
'''
example for using the rem3d api client: fetch model evaluations and cross-sections

python api_model_explorer api_key

api_key is optional if api_initialization has been run.

(this example will include simple plotting for the call results)

'''

from rem3d.api.client import Client as r3d
from rem3d.api.traveltimes import TT
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

    TT_Instance=TT(conn) # give the model instance the connection object

    # list available models & kernels
    model_list=TT_Instance.listModels()
    print("\nModel List")
    print(model_list)


if __name__== "__main__":
    main()

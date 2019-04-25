''' model api class '''
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,requests
import pandas as pd
import numpy as np


class Model(object):

    def __init__(self,client):
        self.r3dc=client
        self.endpoint='/models'
        return

    def listModels(self,args={}):
        '''
        Fetch a list of available models.

        Input parameters:
        ----------------
        none

        Output:
        ------
        result: dictionary containing list of models and model info.

        '''
        args['task']='listModels'
        return json.loads(self.r3dc.call(self.endpoint,dict(args),60))

    def evaluate_points(self,args_in={}):
        '''
        Evaluates list of lat/lon/depth points for a given parameter ('vs')

        Input parameters:
        ----------------
        args_in: dictionary of arguments:

            required args:
                'lat': latitude in degrees
                'lon': longitude in degrees
                'depth': depth in km

                lat, lon and depth are scalars, lists or numpy arrays.
                arrays must be 1d, all must be the same size.

            optional args (default):
                'model': model to use, string ('S362ANI+M')
                'kernel': model to use, string ('BOX25km_PIX1X1')
                'parameter': parameter to fetch, string ('vs')
                'interpolated': 1/0 for interpolation with KdTree (1)

                model+kernel must match a model file.  (see listModels() method)

        Output:
        ------

        result: dictionary with scalar/list/numpy array of same size as input
        '''
        args_in['task']='evaluate_points'
        args=dict(args_in)
        return_numpy=False
        for arg in ['lat','lon','depth']:
            if isinstance(args[arg],np.ndarray):
                return_numpy=True
                args[arg]=args[arg].tolist()
            args[arg]=json.dumps({'vals':args[arg]})
        print(args)

        result=json.loads(self.r3dc.call(self.endpoint,args,5*60))

        if return_numpy:
            param=result['parameter']
            result[param]=np.asarray(result[param])

        return result

    def crossSection(self,args={}):
        '''
        Interpolation of cross-section along a great circle path

        Input parameters:
        ----------------

        args: dictionary of arguments:

            required args:
                'lat': starting latitude in degrees, float,
                'lon': starting longitude in degrees, float,
                'azimuth': direction to move from starting lat/lon
                           (degrees clockwise from North, 90=East, 270=West)
                'gcdelta': great circle distance (degrees) to move along azimuth

            optional args (default):
                'model': model to use, string ('S362ANI+M')
                'kernel': model to use, string ('BOX25km_PIX1X1')
                'parameter': parameter to fetch, string ('vs')
                'interpolated': 1/0 for interpolation with KdTree (1)
                'quickInterp': 1/0 for quick interpolation (0)

            model+kernel must match a model file.  (see listModels() method)

        Output:
        ------

        result: dictionary of numpy arrays

            if args['parameter']='vs',

            result = {'vs': 2d numpy array, 'theta': 1d angle in degrees,
                      'radius': radius in km,
                      'lat': 1d latitude of surface points (not implemented yet),
                      'lon': 1d longitude of surface points (not implemented yet)}
        '''
        args['task']='crossSection'
        json_load = json.loads(self.r3dc.call(self.endpoint,dict(args),5*60))
        parameter=list(json_load.keys())[0]
        results={}
        for keyn in json_load.keys():
            results[keyn]=np.asarray(json_load[keyn])

        return results

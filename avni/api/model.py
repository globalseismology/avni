''' model api class '''
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,requests,os
import pandas as pd
import numpy as np
from configobj import ConfigObj
from .. import tools

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
        args['key']=self.r3dc.key
        args['task']='listModels'
        ModelList=json.loads(self.r3dc.call(self.endpoint,dict(args),60))

        # attach kernel list descriptions
        if ModelList['call_complete']:
            ModelList=self.addConfigDescriptions(ModelList)

        return ModelList

    def addConfigDescriptions(self,ModelList={}):
        '''
        loads kernel, model descriptions from config object
        '''
        filepath = os.path.join(tools.get_configdir(),'attributes.ini')
        if os.path.isfile(filepath):
            # get config file loaded
            parser=ConfigObj(filepath,unrepr=True)

            # compile kernel list to fetch
            kernel_list=[]
            for model in ModelList['3d']['available']:
                kernels=ModelList['3d']['details'][model]['kernel']

                # append kernels to list
                for kern in kernels:
                    if kern not in kernel_list:
                        kernel_list.append(kern)

                # if model in parser['Model3D'].keys():
                #     ModelList['3d']['details'][model]['meta']=parser['Model3D'][model]

            descs={}
            for kern in kernel_list:
                if kern in parser['Kernel_Set'].keys():
                    descs[kern]=parser['Kernel_Set'][kern]

            if len(descs)>0:
                ModelList['kernel_list']=descs

        return ModelList

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
        args_in['key']=self.r3dc.key
        args=dict(args_in)
        return_numpy=False
        for arg in ['lat','lon','depth']:
            if isinstance(args[arg],np.ndarray):
                return_numpy=True
                args[arg]=args[arg].tolist()
            args[arg]=json.dumps({'vals':args[arg]})

        result=json.loads(self.r3dc.call(self.endpoint,args,5*60))
        if result['call_complete']:
            if return_numpy:
                param=result['parameter']
                result[param]=np.asarray(result[param])

        return result

    def depthProfile(self,args_in={},return_numpy=True):
        '''
        Evaluates a depth profile centered on provided lat/lon for a given
        parameter ('vs')

        Input parameters:
        ----------------
        args_in: dictionary of arguments:

            required args:
                'lat': latitude in degrees
                'lon': longitude in degrees

                lat, lon and depth are scalars

            optional args (default):
                'N_depth': depth in km, integer (100)
                'depthMin': min depth, float (0.)
                'depthMax': max depth, float (2890)
                'model': model to use, string ('S362ANI+M')
                'kernel': model to use, string ('BOX25km_PIX1X1')
                'parameter': parameter to fetch, string ('vs')
                'interpolated': 1/0 for interpolation with KdTree (1)

                model+kernel must match a model file.  (see listModels() method)

        Output:
        ------

        result: dictionary with numpy array of parameter, depth
        '''
        args_in['task']='depthProfile'
        args_in['key']=self.r3dc.key
        result=json.loads(self.r3dc.call(self.endpoint,args_in,5*60))

        if result['call_complete']:
            if return_numpy:
                param=result['parameter']
                result[param]=np.asarray(result[param])
                result['depth']=np.asarray(result['depth'])

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
                'includeTopo': 1/0 for incuding topography along transect (0)

            model+kernel must match a model file.  (see listModels() method)

        Output:
        ------

        result: dictionary of results, including numpy arrays

            if args['parameter']='vs',

            result={'parameter':string with parameter name, e.g., 'vs'
                    'vs': 2d numpy array,
                    'depth': 1d array,depth in km,
                    'lat': 1d array, latitude of surface points,
                    'lon': 1d array, longitude of surface points,
                    'theta': 1d array, angular distance along transect [degrees]
                    'topo': 1d array, topography along transect if includeTopo
                    }
        '''
        args['task']='crossSection'
        args['key']=self.r3dc.key
        json_load = json.loads(self.r3dc.call(self.endpoint,dict(args),5*60))

        if json_load['call_complete']:
            results={}
            not_arrays=['parameter','call_complete','topo']
            for keyn in json_load.keys():
                if keyn not in not_arrays:
                    results[keyn]=np.asarray(json_load[keyn])
                else:
                    results[keyn]=json_load[keyn]

            # calculate angular coord (distance along path)
            Nlatlon=results['lat'].size
            results['theta']=np.linspace(0,args['gcdelta'],Nlatlon)

            # adjust topo
            if 'topo' in results.keys() and type(results['topo'])==type(dict()):
                goodTopo=True
                for fld in ['z','lat','lon']:
                    if fld in results['topo'].keys():
                        results['topo'][fld]=np.array(results['topo'][fld])
                    else:
                        goodTopo=False

                if goodTopo:
                    Ntopo=len(results['topo']['z'])
                    results['topo']['theta']=np.linspace(0,args['gcdelta'],Ntopo)

        else:
            results=json_load

        return results

    def fixedDepth(self,args={}):
        '''
        Interpolation at a fixed depth.

        Input parameters:
        ----------------

        args: dictionary of arguments:

            required args:
                'depth': the constant depth to use in km

            optional args:
                'lat1': starting latitude in degrees, float (-90.)
                'lon1': starting longitude in degrees, float (0.)
                'lat2': end latitude in degrees, float (90.)
                'lon2': end longitude in degrees, float (360.)
                if using any of the optional above args, you must use all of
                them. The result will be values in a grid formed by taking
                start and end coordinates are the opposing points of a box
                with the vertex coordinates (lat2,lon1), (lat2,lon2),
                (lat1,lon1), (lat1,lon2). If not specifying coordinates,
                will extract values at a fixed depth for the whole earth.

                'Nlat': number of latitude points to extract, int (100)
                'Nlon': number of longitude points to extract, int (200)
                'model': model to use, string ('S362ANI+M')
                'kernel': model to use, string ('BOX25km_PIX1X1')
                'parameter': parameter to fetch, string ('vs')
                'interpolated': 1/0 for interpolation with KdTree (1)
                'quickInterp': 1/0 for quick interpolation (0)

            model+kernel must match a model file.  (see listModels() method)

        Output:
        ------

        result: dictionary of results, including numpy arrays

            if args['parameter']='vs',

            result={'parameter':string with parameter name, e.g., 'vs'
                    'vs': 2d numpy array with shape (len(lat),len(lon)),
                    'lat': 1d array, latitude
                    'lon': 1d array, longitude
                    }
        '''
        args['task']='fixeddepth'
        args['key']=self.r3dc.key
        json_load = json.loads(self.r3dc.call(self.endpoint,dict(args),5*60))
        results={}

        # convert the results to np arrays
        param=json_load['parameter']
        for fld in [param,'lat','lon']:
            results[fld]=np.array(json_load[fld])

        # copy over some of the other returns
        results['parameter']=param
        results['input_args']=json_load['input_args']

        return results

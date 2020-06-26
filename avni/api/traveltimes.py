''' travel time api class '''
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

class TT(object):

    def __init__(self,client):
        self.r3dc=client
        self.endpoint='/ttimes'
        return


    def listModels(self,args={}):
        '''
        Fetch a list of available models.

        parameters
        -----------
        none

        output
        ------
        ModelList: dictionary containing list of models and model info.

        '''
        args['key']=self.r3dc.key
        args['task']='listModels'
        ModelList=json.loads(self.r3dc.call(self.endpoint,dict(args),60))

        return ModelList

    def predictPaths(self,args={}):
        '''
        result=predictPaths(args)

        predicts travel times for the specified model and paths.

        parameters
        ----------
        args  dictionary, with required and optional keys as follows:

        for model tracing, user must provide one of the following three
        (1) 'start_lat','start_lon','end_lat','end_lon' : start and end
             coordinates in degrees
        (2) 'start_lat','start_lon','azimuth','dist_deg': start coordinates with
            azimuth and great circle distance in degrees.
        (3) 'dist_deg' great circle distance in degrees (FOR 1D tracing ONLY)

        additional required keys:

        'source_depth_km'  source depth in km

        Additional optional parameters

        'type'  tracing type, '1D' or '3D' (default is '1D'), string
        'component'  tracing component, default 'PSV', string
        'model'  the model to use, default is 'PREM_tt_table.h5'. Use
                listModels() to see available models
        'phase' seismic phase to trace, default 'PP'. Use listPhases() to see
                available phases

        output
        ------

        result  dictionary containing results:
        result['ttimes'] travel time in seconds
        '''
        args['key']=self.r3dc.key
        args['task']='predictPaths'

        args = self.processPathArgs(args)

        result=json.loads(self.r3dc.call(self.endpoint,dict(args),60))
        return result

    def processPathArgs(self,args):
        '''
        processPathArgs(self,args)

        processes path arguments in args dict, called by predictPaths()
        '''
        goodArgs=False
        if 'start_lat' in args.keys() and 'end_lat' in args.keys():
            # calculate azimuth, distance
            args['azimuth']=0.
            args['dist_deg']=0.
            goodArgs=True
        elif 'start_lat' in args.keys() and 'azimuth'  in args.keys() and 'dist_deg' in args.keys():
            goodArgs=True
        elif 'dist_deg' in args.keys() and 'start_lat' not in args.keys():
            args['type']='1D'
            args['start_lat']=0.
            args['start_lon']=0.
            args['azimuth']=0.
            goodArgs=True
        else:
            msg=('predictPaths requires: (1) start, end coordiantes or (2) start '
                'coordinates with azimuth and gc distance for 3d models. \n1d '
                'models require (1), (2) or can simply provie gc distance')
            raise ValueError(msg)

        return args

    def listPhases(self,args={}):
        '''
        listPhases()
        lists available phases for component of a given model. Use listModels()
        to see list of available models.

        parameters
        ----------
        args   dictionary with optional keys:
            'type'   '1D' or '3D', default '1D'
            'model'  model to check for phases, default is 'PREM_tt_table.h5'

        output
        ------
        result  dictionary including keys:

        result['component_list'] list of available components
        result['SH'] list of available phases for SH component (if applicable)
        result['PSV'] list of available phases for PSV component (if applicable)        

        '''
        args['key']=self.r3dc.key
        args['task']='listPhases'
        result=json.loads(self.r3dc.call(self.endpoint,dict(args),60))
        return result

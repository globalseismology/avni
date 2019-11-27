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

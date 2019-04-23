''' model api class '''
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,requests
import pandas as pd

class Model(object):

    def __init__(self,client):
        self.r3dc=client
        self.endpoint='/models'
        return

    def listModels(self,args={}):
        args['task']='listModels'
        return json.loads(self.r3dc.call(self.endpoint,args,60))

    def evaluate_at_point(self,args={}):
        args['task']='evaluate_at_point'        
        return json.loads(self.r3dc.call(self.endpoint,args,5*60))

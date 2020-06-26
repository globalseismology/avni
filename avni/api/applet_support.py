''' api support for web_applets. Not recommended for general use. '''
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,requests
import pandas as pd

class App(object):
    ''' generic app support class '''
    def __init__(self,client):
        self.r3dc=client
        return

class SW(App):
    ''' surface wave specific app support '''

    def __init__(self,client):
        App.__init__(self,client)
        self.endpoint='/surfacewaves'
        return

    def buildCommonData(self,args):
        ''' SWbuildCommonData(self,args)
                for SW app. builds events common to two SW data for user selections
                returns filename (file stored in cache)

            args={
                'grp_1': group 1 abbreviation
                'grp_2': group 2 abbrevation
                'period': period to compare at. string, e.g., '100.0'
                'overtone':overtone to compare at. string, e.g., 0
                'wave':wave type 'R' or 'L'
                'orbit': orbit 1, 2, 3. concanenated with wave to get R1, R2, etcself.
                'dat_type': 'processed' or 'raw'
                }

            selections must exist in HDF keys.
        '''
        args['task']='buildCommonData'
        args['key']=self.r3dc.key
        return json.loads(self.r3dc.call(self.endpoint,args,60*5))

    def filterCommonData(self,args):
        ''' SWfilterCommonData(self,args)
            for SW app: args includes plotly/dash graphical selections. See
                        avni-applets.
            returns filename (file stored in cache)
        '''
        args['task']='filterCommonData'
        args['key']=self.r3dc.key
        return json.loads(self.r3dc.call(self.endpoint,args,60*5))

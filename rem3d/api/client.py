''' client for accessing rem3d api '''

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,requests
import pandas as pd

class Client(object):

    def __init__(self,api_key,timeout=300):
        self.timeout=timeout
        self.key=api_key
        self.checkConnection()
        return

    def checkConnection(self):
        ''' checks if the rem3d api is accessible'''
        self.base_url=None
        for url in ['http://127.0.0.1:5000']:

            try:
                r = requests.get(url, timeout=self.timeout)
            except:
                r = None

            if r is not None:
                print("rem3d api is live at ")
                print(url)
                self.base_url=url
                stats=self.checkUserStats()
                print(stats['message'])


        if self.base_url is None:
            print("no connection to rem3d api")
        return

    def call(self,path,parameters=None,time_out=None):
        ''' general form for api calls, returns json as python dictionary '''
        url=self.base_url+path
        if time_out is None:
            time_out=self.timeout
        r = requests.get(url,params=parameters,timeout=time_out)
        return r.json()

    def checkUserStats(self):
        '''
        checks stats for a user
        '''
        args={'key':self.key}
        output=json.loads(self.call('/checkuserstats',args))
        return output

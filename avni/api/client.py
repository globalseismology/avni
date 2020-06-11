''' client for accessing avni api '''

# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,requests
import pandas as pd
from configobj import ConfigObj
from .. import tools
import os

class Client(object):

    def __init__(self,api_key='',timeout=300,api_init_file=None):
        '''
        Client class
        Input:
        ------
            api_key: the api_key to use. If not set, will try to find one
            api_init_file: if api_key is not set, this is the absolute path of
                the api config file. Only needed if setApiConfig() was
                called with api_init_file specified.
            timeout: default timeout in seconds for api call.
        '''
        self.timeout=timeout

        if api_key=='':
            self.key=self.searchForApiKey(api_init_file)
        else:
            self.key=api_key

        self.checkConnection()
        return

    def checkConnection(self):
        ''' checks if the avni api is accessible'''
        self.base_url=None

        connected=False
        for url in ['http://127.0.0.1:5000','http://maurya.umd.edu:41559']:
            if connected is False:
                try:
                    r = requests.get(url, timeout=self.timeout)
                except:
                    r = None

                if r is not None:
                    print("avni api is live at ")
                    print(url)
                    self.base_url=url
                    stats=self.checkUserStats()
                    print(stats['message'])
                    connected=True


        if self.base_url is None:
            print("no connection to avni api")
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

    def searchForApiKey(self,api_init_file=None):
        '''
        searches for API Key
        Input:
        ------
            api_init_file: the api.ini file full path. If None, will default to
                           looking in avni/config/api.ini
        Output:
        -------
            api_key: the api key. Empty string if not found
        '''
        parser=None
        api_key=''
        if api_init_file is None:
            api_init_file = os.path.join(tools.get_configdir(),'api.ini')

        if os.path.isfile(api_init_file):
            parser=ConfigObj(api_init_file)
            api_key=parser['API']['key']

        return api_key

    def setApiConfig(self,api_key='',api_init_file=None):
        '''
        builds the api.ini file
        Input:
        ------
            api_key: the api key to store. Will pull from self.key if not
                     specified.
            api_init_file: the file to store key (subsequents calls to api will
                            need to specify the directory if not using the
                            default). Default is avni/config/api.ini
        Output:
        -------
            no output
        '''

        if api_key=='':
            if self.key=='':
                msg="No api key to set! Initialize Client with api_key or  "
                msg=msg+" call setApiConfig with api_key='APISTRING' "
                print(msg)
            else:
                api_key=self.key

        if api_init_file is None:
            api_init_file = os.path.join(tools.get_configdir(),'api.ini')

        if os.path.isfile(api_init_file):
            print(api_init_file+' already exists. Remove it to reset.')
        else:
            config = ConfigObj()
            config.filename = api_init_file
            config['API'] = {}
            config['API']['key'] = api_key
            config.write()
            print("API key stored in " + api_init_file)

        return

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

class f2pyWrapper(object):
    ''' f2py wrapper '''

    def __init__(self,client):
        self.r3dc=client
        self.endpoint='/f2py'
        return

    def listf2py(self,args={}):
        '''
        lists available f2py functions

        Input:
        -----
        args: None

        Output:
        ------
        result: dictionary containing list of f2py functions
        '''
        args={'key':self.r3dc.key}
        args['f2pyfun']='listfunctions'
        args['f2pyArgs']=json.dumps({})
        f2pyList=json.loads(self.r3dc.call(self.endpoint,args,60))

        return f2pyList

    def callf2py(self,f2pyfun='f2pyfun',f2py_args={}):
        '''
        Abstract function for calling f2py functions in a general way

        Input parameters:
        ----------------
        args: dictionary with required keys:
            f2pyfun: the function to call (string, required)
            f2py_args: the arguments for the f2py function as a kwdict.(required)

        Output:
        ------
        result: result of the f2py function call

        '''
        args={'key':self.r3dc.key}
        args['f2pyfun']=f2pyfun
        args['f2pyArgs']=json.dumps(f2py_args)

        f2pyResult=json.loads(self.r3dc.call(self.endpoint,args,60))

        return f2pyResult

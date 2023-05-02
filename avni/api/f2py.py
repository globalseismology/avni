# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,os
import numpy as np
from configobj import ConfigObj
from .. import tools

class f2pyWrapper(object):
    ''' f2py wrapper '''

    def __init__(self,client):
        self.r3dc=client
        self.endpoint='/f2py'
        self.configfi=os.path.join(tools.get_installdir(),'config','f2py.ini')
        self.configParse=ConfigObj(self.configfi,unrepr=True)
        return

    def listf2py(self):
        '''
        lists available f2py functions in config/f2py.ini

        Input:
        -----
        args: None

        Output:
        ------
        result: dictionary of avaialable f2py functions & arguments
        '''
        f2pyList={}
        print_keys=['description','arguments','output']
        for func in self.configParse.keys():
            f2pyList[func]={}
            for pkey in print_keys:
                f2pyList[func][pkey]=self.configParse[func][pkey]

        return f2pyList

    def callf2py(self,function='f2pyfun',args={}):
        '''
        Abstract function for calling f2py functions in a general way

        Input parameters:
        ----------------
        function: the function to call (string, required)
        args: kwdict for f2py function (required if function has input)

        Output:
        ------
        result: result of the f2py function call

        '''

        if function in self.configParse.keys():
            api_args={'key':self.r3dc.key}
            api_args['f2pyfun']=function
            api_args['f2pyArgs']=self.jsonF2pyArgs(args,function)

            result=json.loads(self.r3dc.call(self.endpoint,api_args,60))
            f2pyResult=self.formatResult(result,function)
        else:
            f2pyResult=function+" is not a valid f2py function"

        return f2pyResult

    def jsonF2pyArgs(self,args,function):
        '''
        converts f2py arguments into json-able form, stores f2py positional order
        and expected type. returns json of the dict.
        '''

        this_func=self.configParse[function]

        args_out={}
        for arg in args.keys():

            # store the argument value, convert types if necessary
            argval=args[arg]
            args_out[arg]={}
            if this_func['arguments'][arg] in ['nparray','nparray_F']:
                args_out[arg]['value']=argval.tolist()
            else:
                args_out[arg]['value']=argval

            # store the required position and expected type of argument
            args_out[arg]['pos']=this_func['argpos'][arg]
            args_out[arg]['type']=this_func['arguments'][arg]


        jsond_args=json.dumps(args_out)
        return jsond_args

    def formatResult(self,result,function):
        ''' returns f2py result in expected tuple, preserving positioning '''
        f2pylist=[]
        this_func=self.configParse[function]

        for iarg in range(0,len(result)):
            # find the name for the output arg at this position
            thisarg=None
            for argkey in this_func['outputtypes'].keys():
                if this_func['outpos'][argkey]==iarg:
                    thisarg=argkey

            if thisarg is None:
                print("could not find output arg")
            else:
                # convert if necessary, append to list
                if this_func['outputtypes'][argkey] in ['nparray']:
                    result[iarg]=np.array(result[iarg])
                f2pylist.append(result[iarg])

        return tuple(f2pylist)

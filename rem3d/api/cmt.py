''' global cmt api class '''
# python 3 compatibility
from __future__ import absolute_import, division, print_function
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import *

# imports for client:
import json,requests
import pandas as pd

class CMT(object):

    def __init__(self,client):
        self.r3dc=client
        self.endpoint='/globalcmt'
        return

    def fetchCMTevents(self,cmt_filters={}):
        ''' fetches events from CMT, returns pandas dataframe.
            Inputs:
                cmt_filters = dictionary for filtering:
                    keys can be: 'date_start','date_end','depth_min_km',
                    'depth_max_km','mag_min' or 'mag_max'
                    date keys are strings: 'YYYY-MM-DD' or 'YYYY-MM-DD HH:MM:SS'
                    other keys are integer or float values.
                    API defaults mag_min to 5, no other defaults.
        '''
        cmt_filters['key']=self.r3dc.key
        ds_json=self.r3dc.call(self.endpoint,cmt_filters,60) # returns column-oriented json string
        events=json.loads(ds_json) # load json as python dict
        if 'call_incomplete' not in events.keys():
            events=pd.DataFrame.from_dict(events).reset_index() # convert to dataframe
        return events

#!/usr/bin/env python
'''
example for using the avni api client: filter the CMT catalogue and plot result

python api_fetchCMT_events api_key

api_key is optional if api_initialization has been run.
'''

from avni.api.client import Client as r3d
from avni.api.cmt import CMT
import matplotlib.pyplot as plt
import argparse

def main():
    # parse
    parser = argparse.ArgumentParser(description='fetches CMT events via API')
    parser.add_argument('-k', '--key', type=str,default='',help='api key')
    arg = parser.parse_args()

    # initialize API connection
    if arg.key=='':
        conn=r3d()
    else:
        conn=r3d(api_key=arg.key) # the connection object

    CMTinstance=CMT(conn)

    filters={'min_depth_km':15}

    print("fetching events")
    events=CMTinstance.fetchCMTevents(filters)

    print("plotting events")
    plt.plot(events['ep.depth'],events['centroid.cmtdepth'],'.k')
    plt.xlabel('reported depth')
    plt.ylabel('centroid depth')
    plt.show()


if __name__== "__main__":
    main()

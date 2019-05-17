#!/usr/bin/env python
'''
example for using the rem3d api client: checks api stats and status for a user

python api_checkUserCallCount api_key

api_key is optional if api_initialization has been run.
'''
from rem3d.api.client import Client as r3d
import argparse

def main():
    # parse
    parser = argparse.ArgumentParser(description='check API user stats')
    parser.add_argument('-k', '--key', type=str,default='',help='api key')
    arg = parser.parse_args()

    # initialize API connection
    if arg.key=='':
        conn=r3d()
    else:
        conn=r3d(api_key=arg.key) # the connection object

    print(conn.checkUserStats()) # print the stats for the user


if __name__== "__main__":
    main()

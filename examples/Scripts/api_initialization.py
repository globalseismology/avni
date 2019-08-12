#!/usr/bin/env python
'''
example for using the rem3d api client: initializing API key storage

python api_initialization api_key

where api_key is the user's api key

stores the api key default api.ini config file, after which subsequent calls to
initialize Client object do not need api_key argument.

'''

from rem3d.api.client import Client as r3d
import argparse

def main():
    # parse
    parser = argparse.ArgumentParser(description='f2py calls via API')
    parser.add_argument('-k', '--key', type=str,default='',help='api key')
    arg = parser.parse_args()

    conn=r3d(api_key=arg.key) # api connection object
    conn.setApiConfig() # store the api key in default config

if __name__== "__main__":
    main()

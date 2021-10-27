#!/usr/bin/env python
'''
example for using the avni api client: initializing API key storage

python api_initialization api_key

where api_key is the user's api key

stores the api key default api.ini config file, after which subsequent calls to
initialize Client object do not need api_key argument.

'''
from avni.api.client import Client as r3d
import argparse

def main():
    # parse
    parser = argparse.ArgumentParser(description='stores api key in config file')
    parser.add_argument('key', type=str,help='api key')
    arg = parser.parse_args()

    conn=r3d(api_key=arg.key) # api connection object
    conn.setApiConfig() # store the api key in default config

if __name__== "__main__":
    main()

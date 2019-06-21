'''
example for using the rem3d api client: initializing API key storage

python api_initialization api_key

where api_key is the user's api key

stores the api key default api.ini config file, after which subsequent calls to
initialize Client object do not need api_key argument.

'''

from rem3d.api.client import Client as r3d
import sys

conn=r3d(api_key=sys.argv[1])
conn.setApiConfig()

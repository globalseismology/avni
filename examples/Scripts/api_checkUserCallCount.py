'''
example for using the rem3d api client: checks api stats and status for a user

python api_checkUserCallCount api_key

where api_key is the user's api key.
'''
from rem3d.api.client import Client as r3d
from rem3d.api.model import Model
import sys

conn=r3d(api_key=sys.argv[1]) # initialize connection
print(conn.checkUserStats()) # print the stats for the user

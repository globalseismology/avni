'''
example for using the rem3d api client: checks api stats and status for a user

python api_checkUserCallCount api_key

api_key is optional if api_initialization has been run.
'''
from rem3d.api.client import Client as r3d
import sys

if len(sys.argv)>1:
    conn=r3d(api_key=sys.argv[1]) # the connection object
else:
    conn=r3d()
print(conn.checkUserStats()) # print the stats for the user

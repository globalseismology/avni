'''
example for using the rem3d api client: filter the CMT catalogue and plot result

python api_fetchCMT_events api_key

where api_key is the user's api key.
'''

from rem3d.api.client import Client as r3d
from rem3d.api.cmt import CMT
import matplotlib.pyplot as plt
import sys


conn=r3d(api_key=sys.argv[1])
CMTinstance=CMT(conn)

filters={'min_depth_km':15}

print("fetching events")
events=CMTinstance.fetchCMTevents(filters)

print("plotting events")
plt.plot(events['ep.depth'],events['centroid.cmtdepth'],'.k')
plt.xlabel('reported depth')
plt.ylabel('centroid depth')
plt.show()

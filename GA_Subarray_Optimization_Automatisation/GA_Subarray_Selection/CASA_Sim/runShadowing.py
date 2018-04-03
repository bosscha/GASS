#!/usr/bin/python

import sys
from arrayConfigurationReport import *

if len(sys.argv) < 2: sys.exit('No cycle specified.')
cycle = sys.argv[3].upper()

if cycle == 'CYCLE5':
    configs = ['ACA-std.cfg', 'C43-1.cfg', 'C43-2.cfg']
else:
    sys.exit('ERROR: CYCLE NOT SUPPORTED')


a=Shadowing()

for config in configs:

    dec1 , shad1 = a.curveShadowingDeclination(config,'1h',-80,+30,3)

    f = open(config.replace('.cfg', '')+'-1h.pickle', 'w')
    data = [dec1,shad1]
    pickle.dump(data,f)
    f.close()

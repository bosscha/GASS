#!/usr/bin/python

"""
Compute and display the LAS  for different array configuration


"""

# sys.path.insert(0,'/home/stephane/git/alma/ALMA/ArrayConfiguration/')

import os

from ArrayConfigurationReportv2 import *



a = LAS()

las1 = a.curveFluxDiskSize('C40-1.cfg','2h', 1. , 30. , 1., -30., cellsize = '0.3arcsec', hourAngle='transit')

f = open('fluxLoss-C40-1.pickle','w')
pickle.dump(las1,f)
f.close()


las2 = a.curveFluxDiskSize('C40-5.cfg','2h', 1. , 30. , 0.5, -30., cellsize = '0.1arcsec', hourAngle='transit')

f = open('fluxLoss-C40-5.pickle','w')
pickle.dump(las2,f)
f.close()




#### Figure
fig = pl.figure()
ax = fig.add_subplot(111)

ax.plot(las1[0],las1[1],'k-')
ax.plot(las2[0],las2[1],'r--')

ax.set_xlim((0., 30))
ax.set_ylim((0.,1.1))

leg = ax.legend(('C40-1','C40-5'), 'upper right')

ax.set_xlabel('Disk size (arcsec)',size= 15)
ax.set_ylabel('Flux (Jy)', size = 15)

## analytical LAS
xx2 = [6.0,6.0]
yy2 = [0.,0.2]

xx1 = [28.9, 28.9]
yy1 = [0., 0.2]

ax.plot(xx1,yy1,'k-')
ax.plot(xx2,yy2,'r--')


pl.savefig("fluxLoss-cy4.png")
pl.show()




#!/usr/bin/python


import UVW as uv



a = uv.UVW("./C43-1-70/sim.tempCasa.ms")
a.plotUVCoverage(-200.,200.,-200.,200.,"Dec = -70", shadow = False, saveFig = True, figfile = 'C43-1-70.png')

a = uv.UVW("./C43-1-23/sim.tempCasa.ms")
a.plotUVCoverage(-200.,200.,-200.,200.,"Dec = -23", shadow = False, saveFig = True, figfile = 'C43-1-23.png')


#!/usr/bin/python


import UVW as uv



a = uv.UVW("./C43-1-23-1m/sim.tempCasa.ms")
a.plotUVCoverage(-200.,200.,-200.,200.,"Dec = -23", shadow = False, saveFig = True, figfile = 'C43-1-23-1m.png')

a = uv.UVW("./C43-1-23-1h/sim.tempCasa.ms")
a.plotUVCoverage(-200.,200.,-200.,200.,"Dec = -23", shadow = False, saveFig = True, figfile = 'C43-1-23-1h.png')
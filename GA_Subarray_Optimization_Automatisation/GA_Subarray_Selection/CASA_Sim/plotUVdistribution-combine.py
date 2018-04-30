#!/usr/bin/python


import UVW as uv



a = uv.UVW("./cy5-1/cy5-1-full.int.ms")
a.plotUVCoverage(-200.,200.,-200.,200.,"Dec = -23", shadow = False, saveFig = True, figfile = 'C43-1-ACA-combine.png')


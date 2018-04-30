#!/usr/bin/python

import pylab as pl
import pickle
import sys
import os


import arrayConfigurationTools as aC


cycle = 'Cycle5'

a = aC.ArrayInfo(cycle+os.path.sep+'ACA-std.cfg')
a.plotAntPos(xmin = -20 , xmax = 20, ymin = -20, ymax = 20, title = '7-m',xtitle = 14, ytitle = 16,figure = '7-m.pdf', showLabels=False)
pl.show()


# a = aC.ArrayInfo('ACA-9-alt.cfg')
# a.plotAntPos(xmin = -20 , xmax = 20, ymin = -20, ymax = 20, title = '7-m-NS',xtitle = 14, ytitle = 16,figure = '7-m--NS.pdf')


a = aC.ArrayInfo(cycle+os.path.sep+'C43-1.cfg')
a.plotAntPos(xmin = -1200 , xmax = 1200, ymin = -1200, ymax = 1200, title = 'C43-1',xtitle = 950, ytitle = 1050,figure = 'C43-1.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-2.cfg')
a.plotAntPos(xmin = -1200 , xmax = 1200, ymin = -1200, ymax = 1200, title = 'C43-2',xtitle = 950, ytitle = 1050,figure = 'C43-2.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-3.cfg')
a.plotAntPos(xmin = -1200 , xmax = 1200, ymin = -1200, ymax = 1200, title = 'C43-3',xtitle = 950, ytitle = 1050,figure = 'C43-3.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-4.cfg')
a.plotAntPos(xmin = -1200 , xmax = 1200, ymin = -1200, ymax = 1200, title = 'C43-4',xtitle = 950, ytitle = 1050,figure = 'C43-4.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-5.cfg')
a.plotAntPos(xmin = -1200 , xmax = 1200, ymin = -1200, ymax = 1200, title = 'C43-5',xtitle = 950, ytitle = 1050,figure = 'C43-5.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-6.cfg')
a.plotAntPos(xmin = -1200 , xmax = 1200, ymin = -1200, ymax = 1200, title = 'C43-6',xtitle = 950, ytitle = 1050,figure = 'C43-6.pdf', showLabels=False)



a = aC.ArrayInfo(cycle+os.path.sep+'C43-7.cfg')
a.plotAntPos(xmin = -6000 , xmax = 6000, ymin = -6000, ymax = 6000, title = 'C43-7',xtitle = 4800, ytitle = 5300,figure = 'C43-7.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-8.cfg')
a.plotAntPos(xmin = -6000 , xmax = 6000, ymin = -6000, ymax = 6000, title = 'C43-8',xtitle = 4800, ytitle = 5300,figure = 'C43-8.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-9.cfg')
a.plotAntPos(xmin = -7000 , xmax = 7000, ymin = -7000, ymax = 7000, title = 'C43-9',xtitle = 5600, ytitle = 6200,figure = 'C43-9.pdf', showLabels=False)

a = aC.ArrayInfo(cycle+os.path.sep+'C43-10.cfg')
a.plotAntPos(xmin = -7000 , xmax = 7000, ymin = -7000, ymax = 7000, title = 'C43-10',xtitle = 5600, ytitle = 6200,figure = 'C43-10.pdf', showLabels=False)

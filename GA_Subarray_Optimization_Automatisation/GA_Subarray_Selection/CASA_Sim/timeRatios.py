import os
import re
import sys
import glob
import math
import numpy
import itertools

cycle = 5
configs = '*'

def readConfigurationFiles():

    if cycle == '' or configs == '': return {}

    cfgFiles = glob.glob(os.path.expanduser('~/AIV/science/ArrayConfiguration/Tools/Cycle'+str(cycle)+'/'+configs+'.cfg'))

    cfgInfo = {}

    for i in range(len(cfgFiles)):

        cfgName = os.path.basename(cfgFiles[i]).replace('.cfg', '')

        f = open(cfgFiles[i])
        fc = f.readlines()
        f.close()

        cfgInfo[cfgName] = {}

        for j in range(len(fc)):
    
            fc1 = fc[j].strip()

            if fc1[0] == '#': continue
            fc1 = fc1.split()

            coordX = fc1[0]
            coordY = fc1[1]
            coordZ = fc1[2]
            padName = fc1[4]

            cfgInfo[cfgName][padName] = {}
            cfgInfo[cfgName][padName]['coord'] = [coordX, coordY, coordZ]

    return cfgInfo

def getBaselineLengths():

    cfgInfo = readConfigurationFiles()
    
    blInfo = {}
    
    for i in cfgInfo.keys():

        blInfo[i] = {}

        for j in itertools.combinations(cfgInfo[i].keys(), 2):

            dist2 = 0
            for k in range(2):
                dist2 += (float(cfgInfo[i][j[0]]['coord'][k]) - float(cfgInfo[i][j[1]]['coord'][k]))**2

            blInfo[i][j] = math.sqrt(dist2)

    return blInfo

def getBaselineStats():

    blInfo = getBaselineLengths()

    for i in sorted(blInfo.keys()):
    
        print i, round(numpy.min(blInfo[i].values()), 1), round(numpy.mean(blInfo[i].values()), 1), round(numpy.max(blInfo[i].values()), 1)

def getConfigurationOverlap():

    blInfo = getBaselineLengths()
    
    for i in itertools.combinations(blInfo.keys(), 2):

        minBL = []
        maxBL = []

        for j in range(2):
            minBL.append(numpy.min(blInfo[i[j]].values()))
            maxBL.append(numpy.max(blInfo[i[j]].values()))

        if maxBL[1] > maxBL[0]:
            cfgComp = 0
            cfgExt = 1
        else:
            cfgComp = 1
            cfgExt = 0

        if re.search('^ACA*', i[cfgComp], re.IGNORECASE) != None:
            diaComp = 7
        else:
            diaComp = 12

        if re.search('^ACA*', i[cfgExt], re.IGNORECASE) != None:
            diaExt = 7
        else:
            diaExt = 12

        overlapComp = len(numpy.where(numpy.array(blInfo[i[cfgComp]].values()) >= minBL[cfgExt])[0])
        overlapExt = len(numpy.where(numpy.array(blInfo[i[cfgExt]].values()) <= maxBL[cfgComp])[0])
        try:
            timeRatio = ( 1.*overlapExt / overlapComp ) * ( 1.*diaExt / diaComp )**2
        except:
            timeRatio = 0

#         print 'compact = '+i[cfgComp]+' extended = '+i[cfgExt]
#         print 'overlap min/max [m]:', round(minBL[cfgExt], 1), round(maxBL[cfgComp], 1)
#         print 'number & fraction overlap [compact]:', overlapComp, round(1.0*overlapComp/len(blInfo[i[cfgComp]].values()), 4)
#         print 'number & fraction overlap [extended]:', overlapExt, round(1.0*overlapExt/len(blInfo[i[cfgExt]].values()), 4)
#         print 'time ratio', round(timeRatio, 3)
#         print ''

        print i[cfgComp], i[cfgExt], round(timeRatio, 2)

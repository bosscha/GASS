#!/usr/bin/python
# -*- coding: latin-1 -*-

"""
Script to be run from arrayCondifuration.py in CASA to produce the simulations and associated results
It read the simulation parameters in param.sim

HISTORY:
    2013.04.16:
        - first shot
        
        
    2013.04.19:
        - add simulation
    
    2013.04.20:
        - write the report
        - add cleaning simulation for the beam resolution
        
    2013.05.09:
        - add comparison with standard array (Cycle 1)
        - add metric1 to compare two uv-distribution.
        - modify the report

    2013.06.23:
        - add a variable (StandardDirectory) for the standard configurations.
        
        
    2013.07.02:
        - add comparison with Cycle 2 including ACA (test)
        
    2013.07.03:
        - read the disk parameter for the component size
        
        
    2013.09.03:
        - adding sidelobe level
        
    2013.09.05:
        - adding spatial resolution (Geometric)
        
        
    2013.09.12:
        - adding weighting for the imaging and updating the report.
        
    2014.01.02:
        - a web method to deal with the -p option (online = yes)
        
        
    2014.02.10:
        - adapt the simanalyze to different spatial resolution to avoid very large image for very fine resolution
        
    2014.02.11:
        - adapt the simanalyze to different spatial resolution to avoid very large image for very fine resolution
        - tested with CASA 4.2.0
        
    2014.02.12:
        - add a check for the image/psf files

    2014.02.13:
        - update simulation method
        
    2014.02.14:
        - fix a bug with the disk size (increase the digit for small disks)
        
    2014.05.28:
        - change the standard configuration to Cycle 2.
        
    2014.05.29
        - add the uv overlapping btw two configurations
        
    2014.06.03:
        - small fix for the format
        
    
    2014.09.30:
        - Read the cycle option
        
    2014.10.06:
        - fix the cycle / standard comparison
        
    2014.10.08:
        - Add Metric to the array configuration (Voronoi tesselation)
        
    2014.10.14:
        - minor update for the report.
        
    2014.10.22:
        - update of voronoiMetric to avoid error if no step Radius is found
        
    2014.11.03:
        - Error check if the periPart is empty in voronoiMetric
        
    2014.11.10:
        - change cycle 3 for testing 
        
    2014.11.14:
        - put official name for Cycle 3
        - fix some error check in voronoiMetric
          
    2014.11.24:
        - add the option of timeRatio to compute the time ratio if combined with a standard array

    2014.11.25:
        - timeRatio computation
        
    2014.12.23:
        - normalization of the uv density to the total number of BL and modifying/adding plots
        

    2014.12.25:
        - add the plot of the histogram of the BL
               
    2014.12.26:
        - creating metric2 for testing a second metric to compare two UV distribution
        - update the histogram BL.
               
    2014.12.27:
        - coding metric2
        
    2014.12.28:
        - updating metric2 with snrUV and snrTot
        
    2015.03.22:
        - adding the computation of the most critical pads using the Voronoi area summed over all the BL for one pad.

    2015.07.22:
        - stable version : we set the Voronoi metrics.
        
    2015.09.24:
        - adding L20, L80 percentile for the configuration.
        
    2015.09.28:
        - add option readParam in the class __init__ to allow not reading for other use.
        - add the casac import to use ia tool in script
        
    2015.10.28:
        - add the Cycle 4 standard configurations
        
    2016.01.29:
        - Adding the AR, LAS estimated in the report from the percentil formulas
        
    2016.06.13:
        - minor update in the header for explanation.
        
        
    2016.08.04:
        - add jason capabilities for the report.
        
    2016.08.05:
        - add Cycle 5 configurations for matching.
        
        
    2016.08.11
        - modify the estimation of the cellsize for standard configuration simulations
        
    
    2016.12.01
        -adding the file with the list of critical antennas and their criticality
        
        
    2017.05.11:
        - adding a -i option for the param file in the call to the script
        
    2017.05.12:
        - updating the naming scheme
        
    2017.05.14:
        - fixing a bug with project name
        
        
    2017.05.16:
        - json version added.
        
        
    2017.07.21:
        - adapting to CASA version >= 5.0.0
    
        
    2017.08.10:
        - fixing bug
    
        
    2018.03.28:
        - path issue
        
RUN:
CASA> sys.path.insert(0,'/home/stephane/git/ALMA/ALMA/ArrayConfiguration/')

casapy -c arrayConfigurationSim.py

Python wrapper:
python arrayConfiguration.py [Options]

=======

"""


__version__="1.5.2@2018.03.28"
__author__ ="ALMA: SL"


import sys
sys.path.insert(0,'./')

import os, os.path
import shutil
import json
import operator
import UVW
import numpy as np
from math import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import pylab as pl
from taskinit import *

RAD2ARCSEC = 3600. * 180. / pi
SPEEDOFLIGHT = 299792458.0

standardDirectory = './'


ia = casac.image()
ms = casac.ms()

class arrayConfigurationSim:
    "Class to run simulations and analysis from a given array configuration "
    
    
    def __init__(self, simulParamFile = 'param.sim', readParam = True):
        
        self.simParamFile = simulParamFile
        
        if readParam:
            
            self.readParam()
            self.StandardArray = []
            self.setStandardArray()
            self.readPadName()
        
    
    
    def readParam(self):
        "Read the parameters in the param.in file"
        
        f = open(self.simParamFile, "r")
          
        for data in f:
            dataSpl = data.split()
            
            if dataSpl[0] == 'file':
                self.fileConfig = dataSpl[2]
                
            if dataSpl[0] == 'inputfile':
                self.inputfile = dataSpl[2]
                
            elif dataSpl[0] == 'reportFile':
                self.reportFile = dataSpl[2]
                
            elif dataSpl[0] == 'Nant':
                self.numberAntenna = int(dataSpl[2])
                           
            elif dataSpl[0] == 'frequency':
                self.freq = float(dataSpl[2])
                self.frequency = "%sGHz"%(dataSpl[2])
                
            elif dataSpl[0] == 'declination':
                self.declination = float(dataSpl[2])
                        
            elif dataSpl[0] == 'hourAngle':
                self.hourAngle = dataSpl[2]
                
            elif dataSpl[0] == 'duration':
                self.duration = dataSpl[2]
                
            elif dataSpl[0] == 'minBl':
                self.minBl = float(dataSpl[2])
                
            elif dataSpl[0] == 'maxBl':
                self.maxBl = float(dataSpl[2])
                
            elif dataSpl[0] == 'rms':
                self.rms = float(dataSpl[2])
                
            elif dataSpl[0] == 'standard':
                self.standard = dataSpl[2]
                
            elif dataSpl[0] == 'timeRatio':
                self.timeRatio = dataSpl[2]
                
            elif dataSpl[0] == 'disk':
                self.disk = float(dataSpl[2])
            
            elif dataSpl[0] == 'weighting':
                self.weight = dataSpl[2]
                
            elif dataSpl[0] == 'cycle':
                self.cycle = int(dataSpl[2])
                
            elif dataSpl[0] == 'online':
                self.online = dataSpl[2]
                
            elif dataSpl[0] == 'criticalPads':
                self.numberCriticalPads = int(dataSpl[2])
                
            elif dataSpl[0] == 'jasonFile':
                self.jasonFile = dataSpl[2]
                
            elif dataSpl[0] == 'projectName':
                self.projectName = dataSpl[2]
                
            elif dataSpl[0] == 'cleanDir':
                self.cleanDir = dataSpl[2]
                                              
        f.close()
    
    
    def setStandardArray(self):
        "Set the standard array according to the cycle"
        
        if self.cycle == 1:
            self.StandardArray=['C32-1.cfg','C32-2.cfg','C32-3.cfg','C32-4.cfg','C32-5.cfg','C32-6.cfg']
            
        elif self.cycle == 2:
            self.StandardArray=['ACA-std.cfg', 'C34-1.cfg','C34-2.cfg','C34-3.cfg','C34-4.cfg','C34-5.cfg','C34-6.cfg','C34-7.cfg']
            
        elif self.cycle == 3:
            self.StandardArray=['ACA-std10.cfg','C36-1.cfg','C36-2.cfg','C36-3.cfg','C36-4.cfg','C36-5.cfg','C36-6.cfg','C36-7.cfg','C36-8.cfg']

        elif self.cycle == 4:
            self.StandardArray=['ACA-std.cfg','C40-1.cfg','C40-2.cfg','C40-3.cfg','C40-4.cfg','C40-5.cfg','C40-6.cfg','C40-7.cfg','C40-8.cfg','C40-9.cfg']

        elif self.cycle == 5:
            self.StandardArray=['ACA-std.cfg','C43-1.cfg','C43-2.cfg','C43-3.cfg','C43-4.cfg','C43-5.cfg','C43-6.cfg','C43-7.cfg','C43-8.cfg','C43-9.cfg','C43-10.cfg']


    def readPadName(self):
        "read the pad names"
        
        self.padName = []
        
        fconf = open(self.fileConfig)
        
        for line in fconf:
              if line[0] != "#":
                  s = line.split()
                  self.padName.append(s[4])
                  
        fconf.close() 
        

    def simulation(self,fileConfig, projectName = 'sim'):
        "Run the actual simulation"
        
        
        if self.declination <= -10:
            decString="J2000 0h00m00s %3dd00m00s"%(self.declination)
        elif self.declination <0 and self.declination>-10:
            decString="J2000 0h00m00s -0%1dd00m00s"%(-self.declination)
        elif self.declination>=0 and self.declination < 10 :
            decString="J2000 0h00m00s +0%1dd00m00s"%(self.declination)
        elif self.declination >=10:
            decString="J2000 0h00m00s +%2dd00m00s"%(self.declination)
                 
        os.system("rm -Rf %s"%(projectName))
        os.system("rm -Rf %s.cl"%(projectName))
        
        diskSize = "%6.4farcsec"%(self.disk)
        
        cl.addcomponent(dir = decString, flux = 1, freq = self.frequency ,shape = "disk", majoraxis = diskSize, minoraxis = diskSize, positionangle = "0deg")
        cl.rename(projectName+".cl")
        cl.done()         
        
        
        ## map size    
        mapS0 = self.disk * 4.
        mapS = "%5.2farcsec"%(mapS0)
        
        default(simobserve)        
        simobserve( 
                   project = projectName, 
                   complist = projectName+".cl", 
                   compwidth = '8GHz', 
                   antennalist = fileConfig, 
                   hourangle = self.hourAngle, 
                   thermalnoise = "",
                   mapsize =  mapS,
                   totaltime = self.duration )
        
        estimatedResolution = 61800. / (self.maxBl * self.freq)      
        
        if self.standard == "yes":
            fileCfgtmp = fileConfig
            
            nhead = fileCfgtmp.find('./')
            
            if nhead > -1 :
                fileCfgtmp = fileCfgtmp[nhead+2:].split('.')[0]
            else :
                fileCfgtmp = fileCfgtmp.split('.')[0]

             
            print  fileCfgtmp
            msName = "%s/%s.%s.ms"%(projectName,projectName,fileCfgtmp)
            print msName
            
            uv = UVW.UVW(msName)
            ruv = uv.distUV()
            maxBl = max(ruv)
            estimatedResolution = 61800. / (maxBl * self.freq)   
            

        Npixel = int(4. * 5. *  self.disk / estimatedResolution)  
        if estimatedResolution > self.disk:
            Npixel =  int(4. * 5. / estimatedResolution)  
            
        cellSize = "%5.4farcsec"%(estimatedResolution / 5.)                
        
        print("estimatedResolution = %f"%(estimatedResolution))
        print("Npixel = %d"%(Npixel))
        print("Cell Size= %s"%(cellSize))
        
        simanalyze(
                   project = projectName,
                   image = True , 
                   imsize = Npixel ,
                   cell = cellSize ,
                   weighting = self.weight,
                   niter = 0)
        
        
    def getMSName(self):
        "Return the MS name of the project"
        
        resSplit = self.fileConfig.split('/')[-1].split(".")   
        msname = '%s/%s.%s.ms'%(self.projectName,self.projectName,resSplit[0])
        return(msname)
    
    
    def getImageName(self):
        "return the image name"
        
        resSplit = self.fileConfig.split('/')[-1].split(".")
        imagename = '%s/%s.%s.image'%(self.projectName,self.projectName,resSplit[0])
        return(imagename)        

    def getNoisyImageName(self):
        "return the noisy image name"
        
        resSplit = self.fileConfig.split('/')[-1].split(".") 
        imagename = '%s/%s.%s.noisy.image'%(self.projectName,self.projectName,resSplit[0])
        return(imagename)        
                        
    def getBeamName(self):
        "return the beam name"
        
        resSplit = self.fileConfig.split('/')[-1].split(".") 
        beamname = '%s/%s.%s.psf'%(self.projectName,self.projectName,resSplit[0])
        return(beamname)        

    def getNoisyBeamName(self):
        "return the noisy beam name"
        
        resSplit = self.fileConfig.split('/')[-1].split(".")
        imagename = '%s/sim.%s.noisy.psf'%(self.projectName,self.projectName,resSplit[0])
        return(beamname) 
    
       
    def analysis(self):
        "To analyze the simulations and  write the report"
        
        ## counting the antennas
        ant = []
        fconf = open(self.fileConfig,"r")
        
        for line in fconf:
              if line[0] != "#":
                  s = line.split()
                  ant.append([float(s[0]), float(s[1])])
                  
        fconf.close()
        
        self.numberAntenna = len(ant)
        
        #####
        ms = self.getMSName()
        
        print("ms .. %s"%(ms))
        
        uv = UVW.UVW(ms)
        ruv = uv.distUV()
        rds = np.sort(ruv)
        
        self.minBlProj = min(ruv)
        self.maxBlProj = max(ruv)
        
        nuv = len(rds)
        i05 = (int)(0.05 * nuv)
        i10 = (int)(0.10 * nuv)
        i20 = (int)(0.20 * nuv)
        i80 = (int)(0.80 * nuv)
        i90 = (int)(0.90 * nuv)
        
        self.L05Proj = rds[i05]
        self.L05Proj2 = np.percentile(rds, 5)
        self.L10Proj = rds[i10]
        self.L20Proj = rds[i20]
        self.L80Proj = rds[i80]
        self.L90Proj = rds[i90]
         
        ## estimated LAS and AR
        
        lambdaFreq = SPEEDOFLIGHT / (self.freq * 1e9)
        k1 =  0.573807
        k2 =  0.98329
        
        self.LASestimated =  RAD2ARCSEC * k2 * lambdaFreq / self.L05Proj
        self.ARestimated  =  RAD2ARCSEC * k1 * lambdaFreq / self.L80Proj
        
        ## Shadowing
        i, self.fracShadow = uv.shadowing()
        
        
        ## Spatial resolution & LAS
        
        self.LAS =  37100. / ( self.minBlProj * self.freq) 
        
        imageName = self.getImageName()
        
        if  not os.path.isdir(imageName) :
            imageName = self.getNoisyImageName()
            
            if not os.path.isdir(imageName):
                print("## image file not found in sim...")
                print("## severe problem.")
                exit()
            
        ia.open(imageName)
        self.beam = ia.restoringbeam()
        ia.close()
        
        
        
        ## Sidelobe levels
        
        beamName = self.getBeamName()
        
        if  not os.path.isdir(beamName) :
            beamName = self.getNoisyBeamName()
            
            if not os.path.isdir(beamName):
                print("## beam file not found in sim ...")
                print("## severe problem.")
                exit()
        
        
        radiusBeam , beamLevel = self.beamLevels(beamName)  
        self.sidelobe = 100.*self.findMaxLevelSideLobe(radiusBeam, beamLevel, self.beam['major']['value'])
        
        
        ## UV radial distribution
        
        NBIN = 50
        NBINhist = 30
        MINr = self.minBlProj * 0.90
        MAXr = self.maxBlProj * 1.10
        
        radius , rho = uv.radialDensity(MINr, MAXr, NBIN)
        
        ## Normalization to the number of BL
        nbl = self.numberAntenna * ( self.numberAntenna - 1.) / 2.
        dr = (MAXr - MINr) / NBIN
        IRho = 0.
        
        print "N antenna .. "
        print self.numberAntenna
        
        for i in range(0,len(radius)):
            IRho += rho[i] * 2.* pi * radius[i] * dr
            
        rho = rho * nbl / IRho
        
        ### plot ..
        
        fig = pl.figure() 
        ax = fig.add_subplot('111')
        ax.semilogy(radius,rho,'b-',drawstyle='steps')
        ax.set_xlabel(r"$R_{UV}$ (meter)")
        ax.set_ylabel(r"$\rho_{UV}$")
        
        fig.savefig('%s.radialUVdensityLog.png'%(self.inputfile))
        pl.close(fig)
        
        
        fig = pl.figure()      
        ax = fig.add_subplot('111')    
        ax.plot(radius,rho,'b-',drawstyle='steps')
        ax.set_xlabel(r"$R_{UV}$ (meter)")
        ax.set_ylabel(r"$\rho_{UV}$")
        
        fig.savefig('%s.radialUVdensity.png'%(self.inputfile))
        pl.close(fig)
        
        print("## radial UV density plot produced.")
        
        
        ## plot the histogram of baseline 
        
        radius , hist = self.histogramBaselines(self.fileConfig, MINr, MAXr, NBINhist)
        
        fig = pl.figure()      
        ax = fig.add_subplot('111')    
        ax.plot(radius, hist,'b-',drawstyle='steps')
        ax.set_xlabel(r"$R_{UV}$ (meter)")
        ax.set_ylabel(r"$N_{BL}$")
        
        fig.savefig('%s.histogramBL.png'%(self.inputfile))
        pl.close(fig)
        
        
        ### Metric #####
        ### Voronoi Tesselation
        
        self.VORrStepPeriSigma, self.VORrStepAreaSigma, self.VORindex , self.VORzero = self.voronoiMetric()
        
        

    def standardArray(self):
        "Compare the input configuration with the standard ones. The simulations have been done before"
        
        print("## Comparison with the standard configurations")
        
        ## input array
        ms = self.getMSName()
        uv = UVW.UVW(ms)
        ruv = uv.distUV()
        
        self.matching  = {}
        self.uvoverlap = {}
        self.snrArr    = {}
        self.snrUVArr  = {}
        self.snrTotArr = {}
        self.arArr     = {}
        self.lasArr    = {}
        self.timeRatioArray = {}
        
        for array in self.StandardArray:
            name = array.split('.')
            msName = "sim%s/sim%s.%s.ms"%(name[0],name[0],name[0])
            uvtarget = UVW.UVW(msName)
            print("## Compare the distribution with %s"%(array))
            compArr, overlap = self.metric1(uv,uvtarget)
            snrArrFrac , snrUVArr, snrTotArr, arArrFrac , lasArrFrac  = self.metric2(uv,uvtarget)
            
            print("## Matching fraction: %3.1f%%  (uv-overlapping [%3.1f m - %3.1f m])"%(compArr,overlap[0],overlap[1]))
            print("## SNR fraction: %3.2f%%  , AR fraction: %3.2f%% , LAS fraction: %3.2f%%"%(snrArrFrac , arArrFrac, lasArrFrac))
            
            self.matching[array] = compArr
            self.uvoverlap[array] = overlap
            self.snrArr[array] = snrArrFrac
            self.snrUVArr[array] = snrUVArr
            self.snrTotArr[array] = snrTotArr
            self.arArr[array]  = arArrFrac
            self.lasArr[array] = lasArrFrac
            
            
            if self.timeRatio == "yes" :
                print("Compute Time Ratio")
                timeRatioArray = self.computeTimeRatio(uv, uvtarget, overlap)
                print("## Time ratio (apply ACA factor if needed) : %3.1f"%(timeRatioArray))
                self.timeRatioArray[array] = timeRatioArray
            
        print("\n\n")

    def computeTimeRatio(self, uv1, uvtarget, overlap):
        "Compute the time ratio to get the same sensitivity in the overlap uv-range"
        
        N1 = self.countUVrange(uv1, overlap[0], overlap[1])
        N2 = self.countUVrange(uvtarget, overlap[0], overlap[1])
        
        if N2 == 0. :
            tr = 0.
        else :
            tr = N2 / N1 
            
        return(tr)
        
    
    def countUVrange(self, uv, uvmin, uvmax):
        "Count the number of uv-points in the uv range"
        
        if uvmin == 0. and uvmax == 0. : 
            return(0.)
        
        n  = len(uv.u)
        ncount = 0
        
        for i in range(0, n):
            ui = uv.u[i]
            vi = uv.v[i] 
            
            uvdist = sqrt(ui*ui+vi*vi)
            if uvdist >= uvmin and uvdist <= uvmax :
                ncount += 1.
                
        return(ncount)
            
        

    def metric1(self,uv1,uvtarget):
        "compare two uv distributions. The target is the reference"
        
        rtarget = uvtarget.distUV()
        maxBl = max(rtarget)*1.05
        minBl = min(rtarget)*0.95
    
        
        ## An exponential polar grid is taken
        
        Nradius = 30
        Nangle = 15

        
        dr = (log(maxBl)-log(minBl)) / Nradius
        rmin = log(minBl)
        dangle = 2.*pi / Nangle
        
        Ntot = Nradius * Nangle
        
        dens1 = np.zeros([Nradius,Nangle])
        denstarget = np.zeros([Nradius,Nangle])
        
        ## target array uv density
        
        n = len(uvtarget.u)
               
        for i in range(0,n):
            ui = uvtarget.u[i]
            vi = uvtarget.v[i]
            
            indexR = int((log(ui*ui+vi*vi)/2.-rmin)/dr)
            indexAngle = int((atan2(ui,vi)+pi)/dangle)
            
            if indexR >= 0. and indexR < Nradius and indexAngle >=0. and indexAngle < Nangle :
                denstarget[indexR,indexAngle] += 1.
                
        normalisation =  1./len(denstarget.nonzero()[0])
        
        ## input array density
        
        n = len(uv1.u)
        
        for i in range(0,n):
            ui = uv1.u[i]
            vi = uv1.v[i]
            
            indexR = int((log(ui*ui+vi*vi)/2.-rmin)/dr)
            indexAngle = int((atan2(ui,vi)+pi)/dangle)
            
            if indexR >= 0. and indexR < Nradius and indexAngle >=0. and indexAngle < Nangle :
                dens1[indexR,indexAngle] += 1.
        
        
        ## Comparison
        densComp = 0.
        
        for i in range(0,Nradius):
            for j in range(0,Nangle):
                if denstarget[i,j] > 0. and dens1[i,j] < denstarget[i,j]:
                    densComp += normalisation * dens1[i,j] / denstarget[i,j]
                elif denstarget[i,j] > 0. :
                    densComp += normalisation
                   


        ## Add a check to give the overlapping UV-range between both configuration
        ## 
        
        minBlTarget = min(rtarget)
        maxBlTarget = max(rtarget)
        
        ruv1 = uv1.distUV()
        minBluv1 = min(ruv1)
        maxBluv1 = max(ruv1)


        if (minBlTarget > maxBluv1 or minBluv1 > maxBlTarget) :
            uv_overlap = [0,0]
        else :
            uv_overlap = [max(minBlTarget,minBluv1), min(maxBlTarget,maxBluv1)]
            
            

        return(densComp * 100., uv_overlap)
        


    def metric2(self, uv1, uvtarget):
        """
        Testing second metric for comparing two UV distribution. Namely how uv matches uvtarget
        Sensitivity and spatial resolution are measured separately
        """
        
        rtarget = uvtarget.distUV()
        maxBl = max(rtarget) * 1.05
        minBl = min(rtarget) * 0.95 
        
        ## An exponential polar grid is taken
        
        Nradius = 30
        Nangle = 15

        
        dr = (log(maxBl)-log(minBl)) / Nradius
        rmin = log(minBl)
        dangle = 2.*pi / Nangle
        
        Ntot = Nradius * Nangle
        
        dens1 = np.zeros([Nradius,Nangle])
        denstarget = np.zeros([Nradius,Nangle])
        
        ## target array uv density
        
        n = len(uvtarget.u)
               
        for i in range(0,n):
            ui = uvtarget.u[i]
            vi = uvtarget.v[i]
            
            indexR = int((log(ui*ui+vi*vi)/2.-rmin)/dr)
            indexAngle = int((atan2(ui,vi)+pi)/dangle)
            
            if indexR >= 0. and indexR < Nradius and indexAngle >=0. and indexAngle < Nangle :
                denstarget[indexR,indexAngle] += 1.
        
        ## input array density
        
        n = len(uv1.u)
        
        for i in range(0,n):
            ui = uv1.u[i]
            vi = uv1.v[i]
            
            indexR = int((log(ui*ui+vi*vi)/2.-rmin)/dr)
            indexAngle = int((atan2(ui,vi)+pi)/dangle)
            
            if indexR >= 0. and indexR < Nradius and indexAngle >=0. and indexAngle < Nangle :
                dens1[indexR,indexAngle] += 1.
                      
            
        ### Estimate the sensitivity of the trial configuration in the uv range of the target configuration
        ###
        
        snrFrac = 0

        
        for i in range(Nradius):
            totAngleTarget  = denstarget[i,:].sum()
            totAngle1       = dens1[i,:].sum()
            
            if totAngle1 > totAngleTarget or totAngleTarget == 0.:
                snrFrac += 1.
            else :
                snrFrac += sqrt(totAngle1 / totAngleTarget)
             
        
        snrUV  = sqrt( dens1.sum() / denstarget.sum() )
        snrTot = sqrt ( 1.0 * len(uv1.u) / len(uvtarget.u))
        
        print "total"
        print len(uv1.u)
        print len(uvtarget.u)
        
        
        snrFrac = snrFrac / Nradius
        
            
        
        ## Estimate the spatial resolution improvement
        ##
        
        r1 = uv1.distUV()
        maxBl1 = max(r1) * 1.05
        minBl1 = min(r1) * 0.95 
        
        arFrac = maxBl1 / maxBl   
        lasFrac = minBl / minBl1
        
        return(snrFrac * 100 , snrUV, snrTot, arFrac * 100. , lasFrac * 100.)
        
            
            
        
        

    def beamLevels(self,beamName):
        "Return the maxima levels of the beam vs. radius"
        
        
        ia.open(beamName)
        hdr = ia.summary()
        
        ## Image Size
        Nx=hdr['shape'][0]
        Ny=hdr['shape'][1]
        
        ## Pixel reference
        refx=hdr['refpix'][0]
        refy=hdr['refpix'][1]
        
        ##Increment (radians)
        dx=hdr['incr'][0]
        dy=hdr['incr'][1]
        dx *= RAD2ARCSEC
        dy *= RAD2ARCSEC
        

        ### Maximum levels vs. radius
        NBIN = 30
        RMIN = 0.
        RMAX = (refx-Nx)*dx*0.9
        
        
        dr=(RMAX-RMIN)/NBIN
        radius=np.arange(NBIN)*dr+RMIN
        radMaxLevel=np.zeros(NBIN)
        
        ## loop over the pixels
        
        for i in range(Nx):
            for j in range(Ny):
                xx=(i-refx)*dx
                yy=(j-refy)*dy
                rr=sqrt(xx*xx+yy*yy)
                indexR=floor((rr-RMIN)/dr)
                
                if indexR >=0 and indexR < NBIN:
                    val=ia.pixelvalue([i,j])['value']['value']
                    if abs(val) > abs(radMaxLevel[indexR]):
                        radMaxLevel[indexR]=val
        
        ia.close()
        
        return(radius,radMaxLevel)
    
    
    def findMaxLevelSideLobe(self,radius, radMaxLevel, beamsize):
        "Find the maximum level for the side lobes for a given beamsize"
        
        
        nData = len(radMaxLevel)
        dfdr = np.zeros(nData-1)
        
        for i in range(0,nData-1):
            diff = radMaxLevel[i] - radMaxLevel[i+1]
            if  diff  < 0. :
                dfdr[i] = 0
            else:
                dfdr[i] = 1
                
        
        maxLevel = 0.
        index = -1
        for i in range(0,nData-2): 
            index += 1
            if abs(dfdr[i+1]-dfdr[i]) == 1:
                
                mL = max(abs(radMaxLevel[i+2]),abs(radMaxLevel[i+1]))
                
                if maxLevel < mL and radius[index] > beamsize*0.6:
                    maxLevel = mL
                    
        return(maxLevel)
           
           
    def UVZenithSnapshot(self,fileConfig):
        "compute the UV coverage from the antenna position only"     
        
        
        XY = []
        fconf = open(fileConfig)
        
        for line in fconf:
              if line[0] != "#":
                  s = line.split()
                  XY.append([float(s[0]), float(s[1])])
                  
        fconf.close()
        
        pads = self.padName
        
        
        nant = len(XY)
        self.numberAntenna = nant
    
        ## Compute the UV coverage
        uv = np.zeros((nant*(nant-1),2))
        padsUV = []   ## Pads associated with each UV  points
        
        index = 0
        
        for i in range(nant):
            for j in range(i+1,nant):
                uv[index][0] = XY[i][0] - XY[j][0]
                uv[index][1] = XY[i][1] - XY[j][1]
                padsUV.append([pads[i],pads[j]])
                
                index += 1
        
        for i in range(index,nant*(nant-1)):
            uv[i][0] = -uv[i-index][0]
            uv[i][1] = -uv[i-index][1]
            padsUV.append([padsUV[i-index][0],padsUV[i-index][1]])
                
        
        return(uv, padsUV)
    
            
    def histogramBaselines(self, fileConfig, minR, maxR, NBIN):    
        "Compute the histogram of the BL from the UVZenithSnapshot"  
        
        uv , padsUV = self.UVZenithSnapshot(fileConfig)
        
        rbl = np.zeros(uv.shape[0])
        
        for i in range(uv.shape[0]):
            rbl[i] = sqrt(uv[i][0]*uv[i][0] + uv[i][1]*uv[i][1])
            
        hist , bin_edge = np.histogram( rbl , bins = NBIN, range = (minR, maxR) , density = False)
        
        
        radius = np.zeros(NBIN)
        for i in range(NBIN):
            radius[i] = (bin_edge[i] + bin_edge[i+1]) / 2.
        
        ## The 0.5 factor is to get the total number of BL w/o the symmetrization
        return( radius, hist * 0.5 )
        
    
              
    def voronoi_perimeter(self,ver, region):
        "compute the perimetet  in the voronoi region and vertices"
        
        perimeter = 0.
        sizeRegion = len(region)
        
        
        for index in range(sizeRegion):
            x0 = ver[region[index]][0]
            y0 = ver[region[index]][1]
                
            x1 = ver[region[(index+1) % sizeRegion]][0]
            y1 = ver[region[(index+1) % sizeRegion]][1]
                
            if region[index]!= -1 and region[(index+1) % sizeRegion] != -1 :
                perimeter += sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
                
        return(perimeter)
        
        
    def voronoi_area(self,ver, region):
        "compute the area   in the voronoi region and vertices"
        
        area = 0.
        sizeRegion = len(region)
        
        
        if -1 in region:
            area = 0.
            
        else :
            for index in range(sizeRegion):
                x0 = ver[region[index]][0]
                y0 = ver[region[index]][1]
                
                x1 = ver[region[(index+1) % sizeRegion]][0]
                y1 = ver[region[(index+1) % sizeRegion]][1]
                
                area += 0.5 * (x0*y1-x1*y0)
        
                
        return(abs(area))            


    def criticalPadsVoronoi(self, area , padsUV):
        "Compute the list of the most critical pads by order of priority using the maximum of the area sum for the Voronoi cells"
        
        
        pads = self.padName
        nant = self.numberAntenna
        
        padSum ={}
        for pad in pads:
            padSum[pad] = 0.
    
        index = 0
        for  value in area:
            padSum[padsUV[index][0]] += value
            padSum[padsUV[index][1]] += value
            index += 1
      
        maxSum = max(padSum.iteritems(), key=operator.itemgetter(1))[1]
        
        for k in padSum.iterkeys():
            padSum[k] = 100. * padSum[k] / maxSum        
               
        sorted_padsum = sorted(padSum.items(), key=operator.itemgetter(1), reverse = True) 
        
        return(sorted_padsum)
    

    def voronoiMetric(self):
        "Compute some metrics using Voronoi Tesselation on the array configuration"
        
        print("## Voronoi tesselation ...")
        
        uvZenith , padsUV  = self.UVZenithSnapshot(self.fileConfig)
        ndat = len(uvZenith)
        
        vor = Voronoi(uvZenith)
        
        
        ver = vor.vertices
        reg = vor.regions
        pt  = vor.point_region
        

        radius = np.zeros(ndat)
        peri   = np.zeros(ndat)
        area   = np.zeros(ndat)
        
        #3 Statistics on the regions  wrt. to the radial uv-distance
        
        for i in range(ndat):
            
            xx = uvZenith[i][0]
            yy = uvZenith[i][1]
            
            region = reg[pt[i]]
            
            radius[i] =   sqrt(xx*xx+yy*yy)
            peri[i]   =   self.voronoi_perimeter(ver,region)
            area[i]   =   self.voronoi_area(ver,region)
   
          
        rmax = max(radius)
        print rmax
        
        sigmaPeri = np.zeros(100)
        sigmaArea = np.zeros(100)
        
        ### compute Critical Pads
        
        self.criticalPads = self.criticalPadsVoronoi(area , padsUV)
        
        ### Compute the metrics
        ## radius in percentage of the first  strong step for the std of perimeter
        rStepPeriSigma = 0.
        rStepAreaSigma = 0.
        
        rStepPeriFound = False
        rStepAreaFound = False
        
        ## PArameter to compute the step - threshold change.
        #####
        stepPeriThreshold = 0.2
        stepAreaThreshold = 0.5
        rMinToConsider = 20
    
        rStepPeriSigma = 100.0
        rStepAreaSigma = 100.0
        
        for i in range(1,100):
            ind = np.where(radius < i * rmax / 100.)
            periPart = peri[ind] 
            areaPart = area[ind] 
                        
                
            if len(ind[0]) > 1:            
                sigmaPeri[i] = periPart.std()
                sigmaArea[i] = areaPart.std()
            else :
                sigmaPeri[i] = 0.
                sigmaArea[i] = 0.
            
            if sigmaPeri[i] == 0. and sigmaPeri[i-1] == 0.:
                diffPeriSigma = 0.
            else:            
                diffPeriSigma = 2.0 *  (sigmaPeri[i] - sigmaPeri[i-1]) / (sigmaPeri[i] + sigmaPeri[i-1])
                
            if sigmaArea[i] == 0. and sigmaArea[i-1] == 0.:
                diffAreaSigma = 0.
            else:            
                diffAreaSigma = 2.0 *  (sigmaArea[i] - sigmaArea[i-1]) / (sigmaArea[i] + sigmaArea[i-1])

            if not rStepPeriFound and diffPeriSigma > stepPeriThreshold and i > rMinToConsider:
                rStepPeriFound = True
                rStepPeriSigma = i * 1.0
                
            if not rStepAreaFound and diffAreaSigma > stepAreaThreshold and i > rMinToConsider:
                rStepAreaFound = True
                rStepAreaSigma = i * 1.0
        
        ## Linear fit of Log(Y) and X
        ## We do not use Area because of some 0 value. 
        
        # select Radius
        indPeri = np.where(radius < rStepPeriSigma * rmax / 100. )
        rp1 = radius[indPeri]
        periPart = peri[indPeri]
        
        indPeri2 = np.where(rp1 > rMinToConsider * rmax /100. )       
        rp2 = rp1[indPeri2]
        periPart2 = periPart[indPeri2]
        
        XXp  = rp2 * 100.0 / rmax
        YYp = np.log(periPart2 / rmax) 
        periFit = np.polyfit(XXp, YYp , 1)
        print periFit
        
        
        #########
        #########
        #########
        
        print(" ")
        print rMinToConsider
        print("# stepPeriThreshold = %3.2f"%(stepPeriThreshold))
        print("# stepAreaThreshold = %3.2f"%(stepAreaThreshold)) 
        print("# Radius of the step of sigma perimeter_sum  (%%) : %3.1f"%(rStepPeriSigma))   
        print("# Radius of the step of sigma area  (%%) : %3.1f"%(rStepAreaSigma))   
        print("# Voronoi index  (smaller better, = 100*alpha) : %3.1f"%(100.* periFit[0]))
        print("# Voronoi zero (smaller better) : %4.2f"%(periFit[1])) 
        print("#")
        print("\n") 

        
        # voronoi_plot_2d(vor)
        # pl.show()
        
        return(rStepPeriSigma, rStepAreaSigma, 100.* periFit[0] , periFit[1])
    
    

    def web(self):
        "produce the plots for the web"
        
        if not os.path.exists("./AC") :
            os.system('mkdir AC')
            
        os.system('cd ./AC')
        os.system('rm *.txt')
        os.system('rm *.png')
            
        ## copy report and plots
        os.system('cp ./report/%s ./AC/ac.txt'%(self.reportFile))
        
        uvRadialDensityFile = "%s.radialUVdensity.png"%(self.inputfile)
        os.system("cp ./report/%s ./AC/uvRadialDensity.png"%(uvRadialDensityFile))
        
        observeFile = "%s.observe.png"%(self.inputfile)
        os.system("cp ./report/%s ./AC/observe.png"%(observeFile))
    
    
    def createCriticalPadFile(self, outfile):
        "Create the critical file pads"
        
        fout = open(outfile,'w')
        txtout = ''
        
        for i in range(min(self.numberCriticalPads,self.numberAntenna)):
            txtout += "%s %4.6f  \n"%(self.criticalPads[i][0], self.criticalPads[i][1])
            
        fout.write(txtout)
        fout.close()   


    def report(self):
        "write report and move the plot"

        if not os.path.exists("./report") :
            os.system('mkdir report')
            
        
        # move the simulation results
        
        os.system('mv sim/sim.tempCasa.observe.png report/%s.observe.png'%(self.inputfile))
        os.system('mv sim/sim.tempCasa.image.png report/%s.image.png'%(self.inputfile))
        os.system('mv %s.radialUVdensity.png report/'%(self.inputfile))
        os.system('mv %s.radialUVdensityLog.png report/'%(self.inputfile))
        os.system('mv %s.histogramBL.png report/'%(self.inputfile))
        
        
        ## create the critical pad file
        criticalPadFile = '%s.criticalPads.txt'%(self.inputfile)      
        self.createCriticalPadFile(criticalPadFile)
        os.system('mv %s report/'%(criticalPadFile))
        
        outJason = {}
        
        fjason = open(self.jasonFile,'w')
        fout  = open(self.reportFile,'w')
        
        
        
        textout =  "=== Array Configuration Summary ==== \n\n"
        
        textout += "Input file : %s \n"%(self.inputfile)
        
        outJason['version']   =  __version__
        outJason['inputCfg']  =  self.inputfile
        
        textout += "Nant : %d \n"%(self.numberAntenna)
        outJason['numberOfAntenna'] = self.numberAntenna
        
        textout += "Declination (Degree): %2.2f \n"%(self.declination)
        outJason['declination'] = float(self.declination)
                         
        textout += "Duration : %s \n"%(self.duration)
        outJason['obsersvingTime'] = self.duration 
        
        textout += "Frequency (GHz) : %2.1f \n"%(self.freq) 
        outJason['frequency'] = self.freq 
        
        textout += "Minimum Proj. BL (meters) : %f \n"%(self.minBlProj)
        outJason['minimumBaseline'] =  self.minBlProj
        
        textout += "Maximum Proj. BL (meters) : %f \n"%(self.maxBlProj)
        outJason['maximumBaseline'] =  self.maxBlProj
        
        textout += "L05, L05bis, L10, L20 Proj. BL (meters): %5.2f, %5.2f, %5.2f, %5.2f \n"%(self.L05Proj,self.L05Proj2,self.L10Proj,self.L20Proj)
        textout += "L80, L90 Proj. BL (meters): %5.2f, %5.2f \n"%(self.L80Proj,self.L90Proj)
        outJason['percentileBaseline'] ={'L05':self.L05Proj,'L10':self.L10Proj,'L20':self.L20Proj,'L80':self.L80Proj,'L90':self.L90Proj}
        
        textout += "RMS (meters) : %5.2f \n"%(self.rms)
        outJason['baselineRMS'] = self.rms 
        
        textout += "Shadowing fraction (%%): %2.1f \n"%(self.fracShadow)
        outJason['shadowing'] = self.fracShadow
        
        textout += "minor beam (arcsec) : %2.4f \n"%(self.beam['minor']['value'])
        textout += "major beam (arcsec) : %2.4f \n"%(self.beam['major']['value'])
        textout += "Spatial resolution (arcsec) : %2.4f \n"%(sqrt(self.beam['minor']['value']*self.beam['major']['value']))
        outJason['spatialResolution']  = {'minor':self.beam['minor']['value'],'major:': self.beam['major']['value'], \
                                               'PA':self.beam['positionangle']['value'],'geometricResolution':sqrt(self.beam['minor']['value']*self.beam['major']['value'])}
        
        textout += "Weighting : %s \n"%(self.weight)
        outJason['weighting'] = self.weight
        
        textout += "Sidelobe level (%%) : %3.1f \n"%(self.sidelobe)
        outJason['sidelobeLevel'] = self.sidelobe
        
        textout += "PA beam (degree): %2.1f \n"%(self.beam['positionangle']['value'])
        textout += "LAS (arcsec) : %2.4f \n"%(self.LAS)
        outJason['LAS'] = self.LAS
        
        textout += "Analytic AR (arcsec) : %2.4f \n"%(self.ARestimated)
        textout += "Analytic LAS (arcsec) : %2.4f \n"%(self.LASestimated)
        textout += "\n"
        
        textout +="Critical Pads [%priority] : \n"
        
        critPads = []
        
        for i in range(min(self.numberCriticalPads,self.numberAntenna)):
            textout += "%s [%2.2f %%]  "%(self.criticalPads[i][0],self.criticalPads[i][1])
            critPads.append([self.criticalPads[i][0],self.criticalPads[i][1]])
            
        outJason['criticalPads'] = critPads
            
            
        textout += "\n\n"
        
        textout += "=== metrics ===\n\n"
        textout += "# Voronoi Perimeter Critical radius (%%) : %2.1f \n"%(self.VORrStepPeriSigma)
        textout += "# Voronoi Area Critical radius (%%) : %2.1f \n"%(self.VORrStepAreaSigma)
        textout += "# Voronoi Compactness  : %4.2f (I: %4.2f, Z: %4.2f) [-1:ok - 0:sparse]\n"%((self.VORindex+self.VORzero)* self.VORrStepAreaSigma / 100., 
                                                                             self.VORindex, self.VORzero)
        
        outJason['Metrics'] = {'voronoiPerimeter':self.VORrStepPeriSigma,'voronoiArea':self.VORrStepAreaSigma, \
                        'voronoiCompactness':(self.VORindex+self.VORzero)* self.VORrStepAreaSigma / 100.}
        
        if self.standard == 'yes':
            textout += "\n"
            textout += "=== Matching with standard configurations ===\n\n"
            
            arrayMatching = []
            
            for array in self.StandardArray:
                textout += "Matching with %12s: %5.1f %%   (uv-overlap: %3.1f m -- %3.1f m)\n"%(array,self.matching[array], self.uvoverlap[array][0],self.uvoverlap[array][1])
                textout += "Matching with %12s: SNR %6.1f %% (uvrange: x%5.2f, total: x%5.2f), AR fraction %6.1f %% , LAS fraction %6.1f %% \n"%(array, self.snrArr[array], self.snrUVArr[array], self.snrTotArr[array], self.arArr[array], self.lasArr[array])
                textout += "\n"
                arrayMatching.append([array,self.matching[array],self.snrArr[array]])
                
                if self.timeRatio == 'yes':
                    textout += "Time ratio with %12s :   %5.2f:1 \n"%(array,self.timeRatioArray[array])
            textout += "(Matching: the higher the better) \n"
            textout += "Note that the SNR is not a direct ratio  but a composite parameter.\n"
            
            outJason['arrayMatching'] = arrayMatching   

        fout.write(textout)
        fout.close()
        

        jas = json.dumps(outJason, indent = 3)
        print >> fjason, jas
        fjason.close()
        
        os.system('mv %s report/'%(self.reportFile))
        
        print(textout)
        
    def cleanDirectories(self):
        "remove some directories to empty disk space"
        
        print("## Cleaning disk space ...")
        
        dirSim = self.projectName
        dirSourceSim = self.projectName+'.cl'
        
        shutil.rmtree(dirSim)
        shutil.rmtree(dirSourceSim)
        
    def run(self):
        "main run method"
        
        self.simulation(self.fileConfig, projectName = self.projectName)
        self.analysis()
             
        
        if self.standard == 'yes':
            
            for array in self.StandardArray:
                name = array.split('.')
                print("## Simulation of the array: %s"%(array))
                self.simulation(standardDirectory+array, projectName = 'sim'+name[0])
            self.standardArray()
        
                
        self.report()
    
        if self.online == "yes" :
            self.web()
    
       
        if self.cleanDir == 'True':
            self.cleanDirectories()
    
    
    
#=============== Main PRogram====================
if __name__ == "__main__":
    
    ## parameter file
    for i in range(len(sys.argv)):
        if sys.argv[i] == '-i':
            paramFile = sys.argv[i+1]
        
    print("## Running simulations ...")  
    print("### Reading simulation parameters in %s"%(paramFile))
    
    a = arrayConfigurationSim(paramFile)
    a.run()
    

    
    
    
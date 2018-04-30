#!/usr/bin/python

"""
Classes to deal with array configuration program acdc for array optimization


INPUT:


HISTORY:
    2012.11.30: 
      - first shot to produce configuration from the best solution
      
    2013.01.11:
      - Adding metrics class to measure quality of an Array
      
    
    2013.01.22:
       - Removed  ArrayConfigurationReportbecause of dependencies of CASA < 3.3 
       - acdc works with CASA 4.0.0
       - simulation rewritten. Only for 12m array sofar.
       - WARNING !! The projectName directory will be removed each time calling simulation !!
       
    2013.01.23:
       - Compute the beam shape for a simulation
      
    2013.01.24:
       - update the createCasaConfiguration to the new acdc output format
       
    2013.01.28:
        - add findMaxLevelSideLobe for the maximum of the sidelobe
        
    2013.01.30:
        - findMaxLevelSideLobe done
        
    2013.02.01:
        - add the quartile metric for the uv-distance distribution
        
    2014.10.01:
        - trim the # comment at the end of the Pad list in the best config
        - run the simulation and some analysis on all .cfg file.
        
    2014.10.02:
        - improving rapport about the simulations.
        - testing a new metric with tesselation (Voronoi or Delaunay) for the uv distribution
        
    2014.10.03:
        - add area to the voronoi tesselation
        - define the metric of the VT
      

    2014.10.06:
        - update of plot for tesselation
      
    
    2014.10.15:
        - add in run the totalTime (mainly for shadowing)
        
    2014.10.22:
        - update of tesselation to avoid error if no step Radius is found
        - completion fraction of the analysis
        
    2014.10.26:
        - minor update
            
    2014.10.30:
        - adding LAS in the report
        
    2014.11.03:
        - changing the size of the image source
        - adding a parameter imsize in run to solve the problem with very large image
        
    2014.11.05:
        - error check for periPart...
        
    2014.11.12:
        - allow voronoi disconnected for ranking
        
    2016.06.22:
        - fix for Casa 4.6
        
        
    2016.07.11:
        - adding messages in analysis
        - add rmaxLevels in beamLevels to speedup computation and restrict to a subimage.
        - check the sidelobe results !!!
      
RUN:

file = 'bestconfigs.txt'
a = acdc(file)
a.run(infoFile = "best-O-1.txt.info", freq = "300GHz")
"""


__version__="0.5.0@2016.07.11"
__author__ ="ALMA: SL"




import arrayConfigurationTools as act
# import ArrayConfigurationReport as acr
import pylab as pl
import os , string
import numpy as np
from casa import simanalyze as sa
from casa import simobserve as so
from casa import componentlist as cl
from casa import image as ia
from casa import table as tb
from math import *
import glob
import UVW

from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay

RAD2ARCSEC=3600.*180./pi


class acdc:
    
    def __init__(self,filename):
        "Filename is the acdc result"
        
        self.filename = filename 
        self.prefix = filename.split('.')[0]
        
        
    def createCasaConfiguration(self):
        "To parse and to create the CASA ALMA configuration file"
        
        
        f = open(self.filename)
        
        indexSolution = 0
        for line in f:
            if line[0] != "#":
                
                ## trim the last comments
                
                ind = line.find("#")
                if ind > 1:
                    line = line[0:ind]
                
                indexSolution += 1
                antenna = line.split()
                
                pads = []
                
                for p in range(4,len(antenna)):
                    pads.append(antenna[p])
                    
                # print pads
                
                name = "%s_%d"%(self.prefix,indexSolution)
                
                aCasa = act.ArrayConfigurationCasaFile()
                aCasa.createCasaConfig(name,listPads = pads)
                
        f.close()
                
                
    def simulation(self,projectName,arrayConfiguration,frequency = '100GHz',totalTime = '30s',dec= -23. , imsizeimage = -1):
        """
        Create the simulations with the following input parameters:
            projectName : name of the project
            arrayConfiguration: name of the CASA file for the array configuration
            frequency : frequency of the observation (Default 100GHz)
            totalTime: total integration time (default 1H)
            dec: declination of the source (default -23 degree)
       
        
        """
        
        print "### Simulation of %s \n"%(arrayConfiguration)
        
        
        if dec <= -10:
            decString = "J2000 0h00m00s %3dd00m00s"%(dec)
        elif dec<0 and dec>-10:
            decString = "J2000 0h00m00s -0%1dd00m00s"%(-dec)
        elif dec>=0 and dec < 10 :
            decString = "J2000 0h00m00s +0%1dd00m00s"%(dec)
        elif dec >=10:
            decString = "J2000 0h00m00s +%2dd00m00s"%(dec)
                 
        os.system("rm -Rf *%s*"%(projectName)) 
        
        cl.addcomponent(dir = decString, flux = 1,freq = frequency,shape="Gaussian",majoraxis = "0.1arcsec",minoraxis = "0.1arcsec",positionangle = "0deg")
        cl.rename(projectName+".cl")
        cl.done()  

            
        so(
                    project = projectName,
                    complist = projectName + ".cl",
                    compwidth = '8GHz',
                    antennalist = arrayConfiguration,
                    totaltime = totalTime,
                    integration = '30s',
                    )
        
        
        sa(
                        project=projectName,
                        image = True , 
                        weighting = "briggs",
                        imsize = imsizeimage
                        # cell = cellsize,
                        # mask = maskClean,
                        # niter = 100,
                        # threshold      =   '0.1mJy',
                        # graphics = 'file'
                        )


        #ia.close()

 
    def beamResolution(self,projectName,arrayConfiguration):
        "Return  the beam resolution (minor,major,PA) of the projectName and arrayConfiguration simulation"
        
        imageName = projectName + "/" + projectName + ".%s.image"%(arrayConfiguration.split('.')[0])      
    
         
        if os.path.exists(imageName):
            ia.open(imageName)
            res = ia.restoringbeam()
            ia.close()
            
            minorAxis = res['minor']['value']
            majorAxis = res['major']['value']
            positionAngle = res['positionangle']['value']
            
            return(minorAxis,majorAxis,positionAngle)
        else:
            return(-1.,-1.,0.)
        
        
    
    def beamLevels(self,projectName,arrayConfigurationPrefix, rmaxLevels):
        "Return the maxima levels of the beam vs. radius"
        
        beamName = projectName + "/" + projectName + ".%s.noisy.psf"%(arrayConfigurationPrefix)
        
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
        RMAX = min((refx-Nx)*dx*0.9, rmaxLevels)
        RMAX2 = RMAX*RMAX
        
        Nx1=int(max(0,refx+floor(rmaxLevels/dx)))
        Nx2=int(min(Nx,refx-floor(rmaxLevels/dx)))    
        Ny1=int(max(0,floor(refy-rmaxLevels/dy)))
        Ny2=int(min(Ny,floor(refy+rmaxLevels/dy)))
        
        print("### Nx1, Nx2, Ny1, Ny2: %d,%d,%d, %d"%(Nx1,Nx2,Ny1,Ny2))
        
        dr=(RMAX-RMIN)/NBIN
        radius=np.arange(NBIN)*dr+RMIN
        radMaxLevel=np.zeros(NBIN)
        
        ## loop over the pixels
        
        for i in range(Nx1,Nx2):
            for j in range(Ny1,Ny2):
                xx=(i-refx)*dx
                yy=(j-refy)*dy
                rr2 = xx*xx+yy*yy
                ## rr=sqrt(xx*xx+yy*yy)
                ## indexR=floor((rr-RMIN)/dr)
                
                if rr2 < RMAX2:
                    rr=sqrt(rr2)
                    indexR=floor((rr-RMIN)/dr)
                    # if indexR >=0 and indexR < NBIN:
                    val=ia.pixelvalue([i,j])['value']['value']
                    if abs(val) > abs(radMaxLevel[indexR]):
                        radMaxLevel[indexR]=val
        
        ia.close()
        
        return(radius,radMaxLevel)
    
    
    def findMaxLevelSideLobe(self,radius, radMaxLevel, beamsize):
        "Find the maximum level for the side lobes"
        
        
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
      
    def arrayCompositePerformance(self, shadowing, sidelobe, minoraxis, majoraxis, vorcriticalradius, vorindex, vorzero):
        "Return a global parameter evaluating the performance of an array configuration"
        
        ranking = 0.
        
        ## shadowwing
        weight1 = 1.
        ranking += (shadowing / 100.) * weight1
        
        ## beam
        weight2 = 2.
        ranking += (majoraxis/minoraxis -1.) * weight2
        
        ## sidelobe
        weight3 = 2.5
        ranking += (sidelobe / 100.) * weight3
        
        ## voronoi tesselation
        weight4 = 1.5
        rankvor = (100. - vorcriticalradius) /100.  + vorindex /20. + pow(10.,vorzero+0.5)
        ranking += rankvor * weight4
        
        ranking /= weight1 + weight2 + weight3 + weight4
        
        finalRanking = 100. - ranking * 100.   ### the higher the  better
        
        return(finalRanking)
                   
                   
                      
        
    def analysis(self, projectName, cfgFilePrefix, frequency, voronoi, rmaxLevels):
        "Analysis and report of the simulation"
        
        report = ""
        report += "# Project : %s \n"%(projectName)
        
        msName = "%s/%s.%s.ms"%(projectName,projectName,cfgFilePrefix)
        imageName = "%s/%s.%s.noisy.image"%(projectName,projectName,cfgFilePrefix)
        
        report += "# MS      : %s \n"%(msName)
        report += "# Image : %s \n"%(imageName)
        
        ## Number of antenna
        print("## info on project..")
        msantenna = msName+"/ANTENNA"
        tb.open(msantenna)
        station = tb.getcol("STATION")
        nant = len(station)
        tb.close
        
        report += "# Number of antennas : %d \n"%(nant)
        report += "# Frequency  : %s \n"%(frequency)
        report += "# total Time : %s \n"%(self.totalTime)
        report += "# declination : %d \n"%(self.declination)
        
        ## UV 
        print("## UV analysis..")
        uv  = UVW.UVW(msName)
        ruv  = uv.distUV()
        r2uv = ruv*ruv
        
        minBlProj = min(ruv)
        maxBlProj = max(ruv)
        rmsBlProj = sqrt( r2uv.sum() / len(ruv))
        
        iGHZ = string.find(frequency,"GHz")
        if iGHZ > 0:
            freq = float(frequency[0:iGHZ])
        else :
            print("## Warning Frequency not found correctly...100 GHz assumed.")
            freq = 100.
        print freq
        
        LASproj =  37100. / ( minBlProj * freq) 
        
        
        report += "# Min Bl (Proj) : %5.1f \n"%(minBlProj)
        report += "# Max Bl (Proj) : %5.1f \n"%(maxBlProj)
        report += "# RMS Bl (Proj) : %5.1f \n"%(rmsBlProj)
        
        
        
        # Shadowing
        print("## Shadowing analysis..")
        i, fracShadow = uv.shadowing()
        
        report += "# Shadowing fraction (%%) : %5.1f \n"%(fracShadow)
        
        # beam resolution
        print("## Beam analysis..")
        
        ia.open(imageName)
        beam = ia.restoringbeam()
        ia.close()
        
        report += "# Minor beam (arcsec) : %2.4f \n"%(beam['minor']['value'])
        report += "# Major beam (arcsec) : %2.4f \n"%(beam['major']['value'])
        report += "# LAS (arcsec) : %5.4f \n"%(LASproj)
        
        
        # sidelobe 
        print("## sidelobe analysis..")
        
        print("### beamLevels..")
        radiusBeam , beamLevel = self.beamLevels(projectName, cfgFilePrefix, rmaxLevels)  
        print("### findMaxLevelSideLobe..")
        sidelobe = 100.*self.findMaxLevelSideLobe(radiusBeam, beamLevel, beam['major']['value'])
        
        report += "# Sidelobe level (%%) : %2.1f \n"%(sidelobe)
        
        
        ### Voronoi Tesselation metric
        print("## Voronoi analysis ..")
        m = metrics()
        if voronoi :
            vor = m.tesselation(msName, showPlot = False)
        else :
            vor = [0.,0.,0.,0.]
        
        report += "# Voronoi Perimeter Critical radius (%%) : %2.1f \n"%(vor[0])
        report += "# Voronoi Area Critical radius (%%) : %2.1f \n"%(vor[1])
        report += "# Voronoi Index  : %4.2f \n"%(vor[2])
        report += "# Voronoi Zero   : %4.2f \n"%(vor[3])
        
        
        ### final ranking
        rank = self.arrayCompositePerformance(fracShadow, sidelobe, beam['minor']['value'], beam['major']['value'], vor[0], vor[2], vor[3] )
        report += "#\n"
        report += "# Final ranking (%% the higher the better)  :  %4.2f \n"%(rank)
        report += "#\n"
        report += "#\n"
        
        
        
        
        return(report, rank)
        
        
    def run(self,infoFile = "none",freq = '100GHz', totalTimeObs = '30s', declination = -23.0, imsize = -1, voronoi = True, rmaxLevels = 1000.):
        "Run the simulations on all the .cfg file"
        
        self.createCasaConfiguration()
        
        listCfg = glob.glob(self.prefix+"*.cfg")
        
        ## info file
        if infoFile  != "none" :
            finfo = open(infoFile,"w")
            
        ## Parameters for simulations
        
        self.totalTime = totalTimeObs
        self.declination = declination
        
        print("## running the simulations")
        print("## patience ...")
        print("##")
        
        bestRank = -1e9
        bestProject = ""
        
        nFile = len(listCfg)
        nAnalyzed = 0.0
        
        for cfgFile in listCfg:
            
            project = cfgFile.split('.')[0]+"_sim"
            
            print("## project %s"%(project))
            
            self.simulation(project, cfgFile, frequency = freq , totalTime = totalTimeObs, dec = declination, imsizeimage = imsize)
            
            print("## analysis...")
            report, rank = self.analysis(project,cfgFile.split('.')[0] , freq, voronoi, rmaxLevels)
            
            if rank > bestRank:
                bestRank = rank
                bestProject = cfgFile
            
            print report
            
            if infoFile != "none":
                finfo.write(report)
               
               
            ## progress
            nAnalyzed += 1.0
            print("### %d files analyzed (%4.1f %%)")%(nAnalyzed, 100. * nAnalyzed / nFile)
            print("#########################")
            
        strOut = "##################"      
        strOut += "## Best Project : \n"
        strOut += "## %s \n"%(bestProject)
        strOut += "## Rank : %4.2f \n"%(bestRank)
        strOut += "## \n\n"
        
        print strOut
          
        if infoFile != "none":
            finfo.write(strOut)
            finfo.close()     


class metrics():
    """
    Several methods to measure metrics on an Array configuration (.cfg)
    """
    
    def __init__(self):
        
        pass
    
        
    def set_Configuration(self,configurationFile):
        "Set configuration file"
        
        self.filename = configurationFile
        self.__readConfiguration__()
        self.nData = len(self.xPos)
        
        
    def __readConfiguration__(self):
        "Read the casa format file"
        
        a = act.ArrayInfo(self.filename)        

        self.xPos = np.zeros(len(a.xPos),float)
        self.yPos = np.zeros(len(a.yPos),float)
        
        for i in range(len(a.xPos)):
            self.xPos[i] = a.xPos[i]
            self.yPos[i] = a.yPos[i]
            
            
    
    def covariance(self):
        "Covariance computation"
        
        covariance = np.cov(self.xPos,self.yPos)
        
        return(covariance)
    
    def dot(self):
        "X,Y dot"
        
        dot = np.dot(self.xPos,self.yPos)
        
        return(dot)

    
    def quartile(self):
        "Return the quartiles (25,50,75) and min,max of the baseline uv-distance distribution" 
        
        a = act.ArrayInfo(self.filename)
        
        print len(a.baseline)
        
        
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
    
        
    def tesselation(self,ms, showPlot = True):
        "Use the tesselation of the UV  points to rank the distribution"
        
        
        uv = UVW.UVW(ms)       
        ndat = len(uv.u)
        dim = (2*ndat,2)   
        points = np.zeros(dim)
        
        print ms
        print dim
        
        for i in range(ndat):
            points[i]  =  [uv.u[i],uv.v[i]]
            points[2*i] = [-uv.u[i],-uv.v[i]]
        
        tri = Delaunay(points)
        vor = Voronoi(points)
        
        ver = vor.vertices
        reg = vor.regions
        pt  = vor.point_region
        

        radius = np.zeros(2*ndat)
        peri   = np.zeros(2*ndat)
        area   = np.zeros(2*ndat)
        
        #3 Statistics on the regions  wrt. to the radial uv-distance
        
        for i in range(2*ndat):
            
            xx = points[i][0]
            yy = points[i][1]
            
            region = reg[pt[i]]
            
            radius[i] =   sqrt(xx*xx+yy*yy)
            peri[i]   =   self.voronoi_perimeter(ver,region)
            area[i]   =   self.voronoi_area(ver,region)

            
          
        rmax = max(radius)
        print rmax
        
        sigmaPeri = np.zeros(100)
        sigmaArea = np.zeros(100)
        
        ### Compute the metrics
        ## radius in percentage of the first  strong step for the std of perimeter
        rStepPeriSigma = 0.
        rStepAreaSigma = 0.
        
        rStepPeriFound = False
        rStepAreaFound = False
        
        ## PArameter to compute the step - threshold change.
        ## 
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
        
        print("# Voronoi tesselation:")
        print rMinToConsider
        print("# stepPeriThreshold = %3.2f"%(stepPeriThreshold))
        print("# stepAreaThreshold = %3.2f"%(stepAreaThreshold)) 
        print("# Radius of the step of sigma perimeter_sum  (%%) : %3.1f"%(rStepPeriSigma))   
        print("# Radius of the step of sigma area  (%%) : %3.1f"%(rStepAreaSigma))   
        print("# Voronoi index  (smaller better, = 100*alpha) : %3.1f"%(100.* periFit[0]))
        print("# Voronoi zero (smaller better) : %4.2f"%(periFit[1])) 
        print("#")
        print("\n") 

        if showPlot:
            # pl.triplot(points[:,0], points[:,1], tri.simplices.copy())
            voronoi_plot_2d(vor)
            
            ##
            fig = pl.figure()
            ax = fig.add_subplot(111)
            ax.plot(range(100),sigmaPeri,marker = "+",linestyle = "-")
            ax.plot(range(100),sigmaArea,marker = "o",linestyle = "--")
            ax.plot(radius * 100.0 / rmax, peri,"g.")
            ax.plot(radius * 100.0 / rmax, area ,"r.")
        
            ax.set_yscale('symlog')
            ax.set_ylim(0.01,max(area))
            pl.show()
            
            print("# Press a key to continue ...")
            a = raw_input()


        return([rStepPeriSigma,rStepAreaSigma, 100.* periFit[0],periFit[1]]  )
        
        
        
        
    
#================================================================================              
if __name__=="__main__":
    " main program"      
    
    file = 'best-O-1.txt'
    a = acdc(file)
    a.run(infoFile = "best-O-1.txt.info", freq = "300GHz")
    
    
  
  
                
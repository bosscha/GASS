#!/usr/bin/python

"""
Classes to create a complete report of  an array configuration with the plots+report.
First it will create the necessary simulations. Then the missing plots+reports


INPUT:

- Antenna List (casa format)
- Band List  : nominal frequencies for each band will be defined
- Type of report: Complete, Normal, Fast  (to be clearly defined)

(Not Yet Implemented)

HISTORY:
    2011.04.03:
        - Definition of the program with the classes.
        
    2011.04.04:
        - Report class
        
    2011.04.17:
        - Add the Shadowing class
        
    2011.04.18:
        - fixing curveShadowingDeclination. Ideally we don't want cleaning but only the UV-coverage.
        - plotShadowingDeclination (+pickle)
        
        
    2011.05.05:
        - Adding the figures class
        
    2011.11.04:
        - adding the class ArrayConfigurationCasaFile
        
    2011.11.06:
        - move the class ArrayConfigurationCasaFile to the file arrayConfigurationFile
        
    2012.03.05:
        - Fix bug in curveShadowDeclination for CASA 3.3
        - Add a method to simulate an image with various components (Gaussian or Disk)
        
    2012.04.04:
        - Add method to analyze the UV-coverage between different array configurations
        
        
    2012.04.11:
        - Add filtering class for the simulation
        
        
    2012.04.19:
        - Modify filtering class to have a larger mapsize and to sue the new task to replace simdata (simobserve , simanalyze)
        
     
    2012.05.17:
        - modify filtering to use a uniform disk rather than a Gaussian
        
        
    2012.05.24:
        - add the beam class to analyze it
        
    2012.08.27:
        - add hourAngle parameter to the Shadowing method.
        
    2012.11.14:
        - Adding the LAS measurement (at 100 GHz) to the curveBeamResolutionDeclination method
        
        
    2013.06.09:
        - adapting arrayConfigurationReport to CASA 4.0.1
        
        
    2013.06.22:
        - Adding LAS class to compute it in function of disk size...
        
    
    2013.07.02:
        - adapt radialDensityArrayConfigurations to CASA 4.0.1
        
    2014.11.07:
        - Adding the script to the AIV CVS
        
    2014.11.08:
        - Important : spatial resolution is computed with Briggs 
        
    2015.07.27:
        - Update the curveFluxDiskSize
        
    2015.10.08:
        - change the LAS definition (MRS) using L05 definition from SCIREQ-328
        IMPORTANT : if MRS_new > MRS_old we keep MRS_old 
        
    2016.03.08:
        - change the mask definition in the LAS class.
        
    2016.03.09:
        - improve LAS class (setting, mapsize, etc)
        
RUN:


CASA> : sys.path.insert(0,'/home/stephane/git/ALMA/ALMA/ArrayConfiguration/')

"""

__authors__="ALMA: SL"
__version__="0.6.2@2016.03.09"


import UVW as uvw
import beamAnalysis as ba
import math, string
import numpy as np
import os
import pickle
from  casa import componentlist as cl
# import simdata as sim
import pylab as pl
from casa import image as ia
from casa import regionmanager as rg
import simanalyze as sa
import simobserve as so
import imstat as ims

import arrayConfigurationTools as aT

RAD2ARCSEC = 206264.9

# execfile('UVW.py')
# execfile('BeamAnalysis.py')

band={'3':'100GHz','6':'230GHz','7':'345GHz','9':'675GHz'}



    
class filtering:
    """
    Simulation of Gaussian Beams (from R. Indebetouw)
    """
    
    def __init__(self):
        
        pass
        
    def runFiltering(self,antennaCfg,trackDuration,frequency,dec,imageTotalSize,resolutionStart,resolutionEnd,resolutionStep):
        """
        
        use the disk100.fits image of a uniform disk of 100 pixel size (diameter)
        
        Run a set of simulations to account for the filtering at different scales.
            antennaCfg: antenna configuration file
            trackDuration : trackDuration
            resolutionStart,resolutionEnd,resolutionStep: range for the resolution to simulate
            dec: declination
            
        OUTPUT:
                resolution, flux: resolution and flux output (1 Jy in entry)
                
                
        """
        
        maskClean = [70,70,180,180]
        
        projectName="tempCurveFiltering"
        
        nStep=math.floor((resolutionEnd-resolutionStart)/resolutionStep)
        resolutionArr=np.arange(nStep)*resolutionStep+resolutionStart
        flux=np.zeros(nStep)
        boxStat = '%d,%d,%d,%d'%(maskClean[0],maskClean[1],maskClean[2],maskClean[3])
        
        os.system("rm -Rf *%s*"%(projectName))
        
        if dec <= -10:
            decString="J2000 10h00m00s %3dd00m00s"%(dec)
        elif dec<0 and dec>-10:
            decString="J2000 10h00m00s -0%1dd00m00s"%(-dec)
        elif dec>=0 and dec < 10 :
            decString="J2000 10h00m00s +0%1dd00m00s"%(dec)
        elif dec >=10:
            decString="J2000 10h00m00s +%2dd00m00s"%(dec)
              
              
            # simulayion with one arcsec component
            ##os.system("rm -Rf *%s*"%(projectName))
            
        print decString
            
        index=0
            
        print resolutionArr
        
        for res in resolutionArr:
                 
            resolutionDisk = "%4.2farcsec"%(res/100.)
            print resolutionDisk
            print antennaCfg
            
            
            so.sim_observe(
                        project = projectName+"%d"%(index),
                        skymodel           =  'disk100.fits',
                        inbright           =  '1Jy/pixel',
                        indirection        =  decString ,
                        incell             =  resolutionDisk,
                        incenter           =  frequency,
                        inwidth            =  '8GHz',
                        antennalist = antennaCfg,
                        totaltime = trackDuration,
                        integration = '10s',
                        direction = decString,
                        maptype = 'ALMA',
                        # mapsize = ['2arcmin','2arcmin'],
                        # pointingspacing = "2arcmin"
                        )
            
            
            sa.sim_analyze(
                        project=projectName+"%d"%(index),
                        image = True , 
                        imsize = imageTotalSize,
                        # cell = cellsize,
                        mask = maskClean,
                        niter = 2000,
                        threshold      =   '0.1mJy',
                        graphics = 'file'
                        )
                       
        
            imageName=projectName+"%d/"%(index)+projectName+"%d.%s.image"%(index,antennaCfg.split('.')[0])
            # imageName=projectName+"%d/"%(index)+projectName+"%d.image"%(index)
            
            print "Read %s"%(imageName)
        
            ia.open(imageName)
            # flux0=ia.statistics()['flux'][0]
            
            stat = ims.imstat(
                               imagename = imageName,
                               box = boxStat ,
                               )
            
            flux0 = stat['flux'][0]
            print flux0
            
            if res == resolutionArr[0]:
                fluxNormalization = flux0
            
            fluxRes = flux0/fluxNormalization
        
            print "Flux: %4.3f"%(fluxRes)
            
            flux[index] = fluxRes
            
            ia.close()
            
            index += 1


        return(resolutionArr,flux)
    



class beam:
    """
    Class to analyze the properties of the beam 
    """

    def __init__(self):
        
        pass
    
    
    def curveBeamResolutiongDeclination(self,antennaCfg,trackDuration,declinationMin,declinationMax,decStep, ha = 'transit'):
        """
        Compute the minor,Major axis vs. Declination for a given configuration
            antennaCfg: antenna configuration
            trackDuration: time duration, e.g. "6h"
            declinationMin,declinationMax: range of declination
            decStep: step for the declination, needs to be an integer...
            
        Output : arrays [declination] [minor] [major]
        """
              
        
        projectName="tempCurveBeamDeclination"
        nStep=math.floor((declinationMax-declinationMin)/decStep)
        decArr=np.arange(nStep)*decStep+declinationMin
        minorAxis = np.zeros(nStep)
        majorAxis = np.zeros(nStep)
        las100 = np.zeros(nStep)
        
        index=0
        
        
        
        # simulation with one arcsec component. We clean the old files.
        os.system("rm -Rf *%s*"%(projectName))
        
        for dec in decArr:
            
            if dec <= -10:
                decString="J2000 0h00m00s %3dd00m00s"%(dec)
            elif dec<0 and dec>-10:
                decString="J2000 0h00m00s -0%1dd00m00s"%(-dec)
            elif dec>=0 and dec < 10 :
                decString="J2000 0h00m00s +0%1dd00m00s"%(dec)
                
            elif dec >=10:
                decString="J2000 0h00m00s +%2dd00m00s"%(dec)
              
              
            # simulayion with one arcsec component
            #os.system("rm -Rf *%s*"%(projectName))
            
            print decString

            #cl.done()
            cl.addcomponent(dir=decString, flux=1,freq='100GHz',shape="Gaussian",majoraxis="0.1arcsec",minoraxis="0.1arcsec",positionangle="0deg")
            cl.rename(projectName+"%d.cl"%(index))
            cl.done()
            
            
            print projectName+"%d.cl"%(index)
            
            so.simobserve(
                    project=projectName+"%d"%(index),
                    complist=projectName+"%d.cl"%(index),
                    compwidth = '8GHz',
                    antennalist=antennaCfg,
                    totaltime=trackDuration,
                    integration = '10s',
                    hourangle = 'ha'
                    )
                       
            
                       
            sa.simanalyze(
                        project=projectName+"%d"%(index),
                        image = True , 
                        weighting = 'briggs',
                        # imsize = 256,
                        # cell = cellsize,
                        # mask = maskClean,
                        niter = 100,
                        # threshold      =   '0.1mJy',
                        # graphics = 'file'
                        )
                       
            psfName= projectName+"%d/"%(index)+projectName+"%d.%s.psf"%(index,antennaCfg.split('.')[0])
            psfQuickName = projectName+"%d/"%(index)+projectName+"%d.%s.quick.psf"%(index,antennaCfg.split('.')[0])
            msName = projectName+"%d/"%(index)+projectName+"%d.%s.ms"%(index,antennaCfg.split('.')[0])
            imageName = projectName+"%d/"%(index)+projectName+"%d.%s.noisy.image"%(index,antennaCfg.split('.')[0])
            
            # if os.path.exists(psfName):
            #     psfRes = psf.psf_calc(psfName)
            # else:
            #     psfRes = psf.psf_calc(psfQuickName)
            
            ia.open(imageName)
            data = ia.restoringbeam()
            ia.close()
                
            print data
            
            minorAxis[index] = data['minor']['value']
            majorAxis[index] = data['major']['value']
            
            print "minor Axis:",minorAxis[index]
            print "major Axis:",majorAxis[index]
            
            ms = uvw.UVW(msName)
            dUV = ms.distUV(noShadow = True)
            duvMin = min(dUV)
            
            rds = np.sort(dUV)
            nuv = len(rds)
            i05 = (int)(0.05 * nuv)
            L05 = rds[i05]
            
            las = RAD2ARCSEC * 0.0017987547479999 / duvMin      ## 0.6 * Lambda / BL_min
            ## Definition of LAS from SCIREQ-328
            las100[index] =  RAD2ARCSEC * 0.983296 * 0.0029979245799998332 / L05
            
#             if las100[index] > las :
#                     las100[index]  = las
#             las100[index]  = las
                    
            print "mininum BL:",duvMin
            print "L05:",L05
            print "LAS (old definition)",las
            print "LAS (100 GHz):",las100[index]
            
            
            
            index+=1
            
        
        # fout=open("curveBeamDec.last","wb")
        # data={'Dec':decArr,'Shadow':fractionShadowing}
        # pickle.dump(data,fout)  
        # fout.close()  
        
        return(decArr,minorAxis,majorAxis,las100)
        





        
class UVcoverage:
    """
    Class to produce the UVcoverage plots and the derived densities, to  provide
    """
    
    def __init__(self):
        
        pass
    
    
    
    def radialDensityArrayConfigurations(self,arrayNameList,totalTimeList,Rmin,Rmax,Nbin,dec):
        """
        Return the radial UV-density for different arrays with a total time of integration for each array 
        INPUT:
        arrayNameList : list of the array to be analyzed.
        totalTimeList:List of total time of integration for each array
        Rmin,Rmax,Nbin: minimum, maximum, bin number for the radius. (common)
        dec: Declination of the source.
        """
    
        
        if dec <= -10:
            decString="J2000 0h00m00s %3dd00m00s"%(dec)
        elif dec<0 and dec>-10:
            decString="J2000 0h00m00s -0%1dd00m00s"%(-dec)
        elif dec>=0 and dec < 10 :
            decString="J2000 0h00m00s +0%1dd00m00s"%(dec)
        elif dec >=10:
            decString="J2000 0h00m00s +%2dd00m00s"%(dec)
                 
        projectName = "tempUVcoverage"
        os.system("rm -Rf *%s*"%(projectName)) 
        
        cl.addcomponent(dir=decString, flux=1.,freq='100GHz',shape="Gaussian",majoraxis="10arcsec",minoraxis="10arcsec",positionangle="0deg")
        cl.rename(projectName+".cl")
        cl.done()         
        
        indexTime = 0
        dens = []
        
        for arrayName in arrayNameList:
            
            print "Configuration: %s \n"%(arrayName)
            
            trackDuration = totalTimeList[indexTime]
            
            so.simobserve(
                    project = projectName,
                    complist = projectName+".cl",
                    compwidth = '8GHz',
                    antennalist = arrayName,
                    totaltime = trackDuration,
                    integration = '10s',
                    )           
            
            msName=projectName+"/"+projectName+".%s.ms"%(arrayName.split('.')[0])
            
            uv = uvw.UVW(msName)
            rad, rhouv = uv.radialDensity(Rmin,Rmax,Nbin)
            
            dens.append(rhouv)
            indexTime += 1
            
        return (rad, dens)
             
    
    
class targetSimulation:
    """
    Class to provide image of a real object seen with the current Array Configuration
    """
    
    def __init__(self):
        
        self.source = 0.
    
    
    def buildImageWithComponentList(self,nameComp,resolutionList,maxSize,xyStep,sh='Gaussian',fluxList=[]):
        """ Input List:
            minRes,maxRes : minimum,maximum of the components
            nComp : number of components 
        """
        
        
        RAString="J2000 0h00m"
        DecString="-30d00m"
        
        xStep = 30.-maxSize/2.
        yStep = 30.+maxSize/2.
        dStep = xyStep
        
        if len(fluxList) != len(resolutionList):
            fluxList=np.zeros(len(resoltuionList),np.int)
            fluxList[:]=1.
        
        indexSource=0
        
        for resolution in resolutionList:
            xStep+=dStep
            if xStep > 30.+maxSize/2.:
                xStep=30.-maxSize/2.
                yStep-=dStep
            fluxSource=fluxList[indexSource]
            indexSource+=1
            
            posString="%s0%3.1fs"%(RAString,xStep/15.)+" %s%3.1fs"%(DecString,yStep)
            print posString, "%5.2farcsec resolution"%(resolution)
            cl.addcomponent(dir=posString, flux=fluxSource,freq='100GHz',shape=sh,majoraxis="%5.2farcsec"%(resolution),minoraxis="%5.2farcsec"%(resolution),positionangle="0deg")
        
        cl.rename(nameComp)
        cl.done()
    


class LAS: 
    """
    Class to study the LAS for a given configuration
    
    """
    
    def __init__(self):

        pass
    
    
    def curveFluxDiskSize(self,antennaCfg, trackDuration, sizeMin, sizeMax, sizeStep, declination, cellsize = '0.2arcsec', 
                          hourAngle='transit', shapeCL = 'Gaussian', threshClean = '0mJy'):
        
        """
        Compute the flux of disk with different size
        
            antennaCfg: antenna configuration
            trackDuration: time duration, e.g. "6h"
            sizeMin,sizeMax: range of the disk size in arcsec. The total flux is 1 Jy
            sizeStep: step for the disk size
            declination: disk declination
            cellsize : size of the cell
        
        Output: array [sizeDisk, flux]
        """
        
        
        ## setting simobserve
        
        projectName="tempDiskFlux"
        nStep = int((sizeMax-sizeMin) / sizeStep)
        sizeArr = np.arange(nStep)*sizeStep + sizeMin
        flux = np.zeros(nStep)
        size = np.zeros(nStep)
        las100 = np.zeros(nStep)
        
        totalSize = 120.
        mapSize = "%3.1farcsec"%(totalSize)
        cellF = float(cellsize[:-6])         ## to be improved ....
        imagesize = int(totalSize/cellF)
        
        index=0
        
        print "image size: %d "%(imagesize)
        print "mapsize :%s"%(mapSize)
        
        # simulation with one arcsec component. We clean the old files.
        os.system("rm -Rf *%s*"%(projectName))
        
        
            
        if declination <= -10:
            decString="J2000 0h00m00s %3dd00m00s"%(declination)
            raDisk  = "0h00m00s"
            decDisk = "%3dd00m00s"%(declination)
            
        elif declination < 0 and dec>-10:
            decString="J2000 0h00m00s -0%1dd00m00s"%(-declination)
            raDisk  = "0h00m00s"
            decDisk = "-0%1d00m00s"%(declination)
            
        elif declination >=0 and dec < 10 :
            decString="J2000 0h00m00s +0%1dd00m00s"%(declination)
            raDisk  = "0h00m00s"
            decDisk = "+0%1d00m00s"%(declination)
            
        elif declination >=10:
            decString="J2000 0h00m00s +%2dd00m00s"%(declination)
            raDisk  = "0h00m00s"
            decDisk = "+0%2d00m00s"%(declination)
              
            
        ## Estimation of the spatial resolution at 100 GHz for the mask
        
        arrCfg = aT.ArrayInfo(antennaCfg)
        
        arrCfg.stats()
        
        maxBL = arrCfg.maxBl
        resEstimated = 61800. / (maxBL * 100.)
        
        print("Estimated RES: %3.2f"%(resEstimated))


        for disksize in sizeArr :
            
            #cl.done()
            cl.addcomponent(dir=decString, flux=1. , freq='100GHz',shape= shapeCL,majoraxis= "%2.2farcsec"%(disksize),minoraxis="%2.2farcsec"%(disksize),positionangle="0deg")
            cl.rename(projectName+"%d.cl"%(index))
            cl.done()
            
            
            print projectName+"%d.cl"%(index)
            
            so.simobserve(
                    project = projectName+"%d"%(index),
                    complist = projectName+"%d.cl"%(index),
                    compwidth = '8GHz',
                    antennalist = antennaCfg,
                    totaltime = trackDuration,
                    integration = '10s',
                    mapsize  = mapSize ,
                    user_pwv = 0. ,
                    )
            
            maxDisk = max(disksize, resEstimated )
                       
            regDiskClean = 'circle[[%s , %s], %3.1farcsec ]'%(raDisk,decDisk, maxDisk * 1.5)
            regDiskFlux  = 'circle[[%s , %s], %3.1farcsec ]'%(raDisk,decDisk, maxDisk * 1.5)
                       

            
            sa.simanalyze(
                        project=projectName+"%d"%(index),
                        image = True , 
                        weighting = 'briggs',
                        imsize = imagesize,
                        cell = cellsize,
                        mask = regDiskClean,
                        niter = 2000,
                        threshold      =   threshClean,
                        imdirection  = decString ,
                        # graphics = 'file'
                        )
                       
            msName = projectName+"%d/"%(index)+projectName+"%d.%s.ms"%(index,antennaCfg.split('.')[0])
            imageName = projectName+"%d/"%(index)+projectName+"%d.%s.noisy.image"%(index,antennaCfg.split('.')[0])

            print(" ")
            print("## analysis ...")
            print "Mask Clean :: %s"%(regDiskClean)
            print "Mask Flux :: %s"%(regDiskFlux)
            print "Component Size :: %2.2f"%(disksize)
            
            ia.open(imageName)
            data = ia.restoringbeam()
            bmaj = data['major']['value']
            
        
            stat = ia.statistics (region = regDiskFlux)
            flux[index] = stat['flux']
            
            ia.close()
                
            print data
            
            ms = uvw.UVW(msName)
            dUV = ms.distUV(noShadow = True)
            duvMin = min(dUV)
            las100[index] = RAD2ARCSEC * 0.0017987547479999 / duvMin      ## 0.6 * Lambda / BL_min
            
            print "mininum BL:",duvMin
            print "LAS (100 GHz & min BL):",las100[index]
            print "Flux: %2.3f"%(flux[index])
            print "--"
            
            
            index+=1
      
        return([sizeArr,flux])
        
        
        
        

class Shadowing:
    """
    Class  for computing shadowing properties
    """
    
    
    def __init__(self):
    
    
        pass
    
    
    
    def curveShadowingDeclination(self,antennaCfg,trackDuration,declinationMin,declinationMax,decStep,hourAngle = 'transit'):
        """
        Compute the shadowing fraction vs. Declination for a given configuration
            antennaCfg: antenna configuration
            trackDuration: time duration, e.g. "6h"
            declinationMin,declinationMax: range of declination
            decStep: step for the declination, needs to be an integer...
            
        Output : array [declination,fraction of shadowing]
        """
 
        projectName="tempCurveShadowingDeclination"
        
        nStep=math.floor((declinationMax-declinationMin)/decStep)
        decArr=np.arange(nStep)*decStep+declinationMin
        fractionShadowing=np.zeros(nStep)
        
        index=0
        
        # simulayion with one arcsec component. We clean the old files.
        os.system("rm -Rf *%s*"%(projectName))
        
        for dec in decArr:
            
            if dec <= -10:
                decString="J2000 0h00m00s %3dd00m00s"%(dec)
            elif dec<0 and dec>-10:
                decString="J2000 0h00m00s -0%1dd00m00s"%(-dec)
            elif dec>=0 and dec < 10 :
                decString="J2000 0h00m00s +0%1dd00m00s"%(dec)
                
            elif dec >=10:
                decString="J2000 0h00m00s +%2dd00m00s"%(dec)
              
              
            # simulayion with one arcsec component
            ##os.system("rm -Rf *%s*"%(projectName))
            
            print decString

            #cl.done()
            cl.addcomponent(dir=decString, flux=1,freq='100GHz',shape="Gaussian",majoraxis="10arcsec",minoraxis="10arcsec",positionangle="0deg")
            cl.rename(projectName+"%d.cl"%(index))
            cl.done()
            
            
            print projectName+"%d.cl"%(index)
            
#             sim.simdata(
            so.simobserve(
                    project=projectName+"%d"%(index),
                    complist=projectName+"%d.cl"%(index),
                    compwidth = '8GHz',
                    antennalist=antennaCfg,
                    totaltime=trackDuration,
                    integration = '10s',
#                     image = False , 
                    hourangle      =  hourAngle,
                    # imsize = 256,
                    # cell = '0.5arcsec',
                    #graphics = 'none'
                    )
                       

            print "Read %s%d.ms"%(projectName,index)
            
            antList=string.split(antennaCfg,'.')
            msName=projectName+"%d/"%(index)+projectName+"%d.%s.ms"%(index,antList[0])
            print msName
            vis=uvw.UVW(msName)
            ii,shad=vis.shadowing()
            fractionShadowing[index]=shad
            
            print "Shadowing fraction: %2.2f"%(shad)
            index+=1
            
        
        fout=open("curveShadowingDec.last","wb")
        data={'Dec':decArr,'Shadow':fractionShadowing}
        pickle.dump(data,fout)  
        fout.close()  
        
        return(decArr,fractionShadowing)
        
    
    def plotShadowingDeclination(self,declination,shadowing):
        "Plot the shadowing fraction vs. the declination"
        
        pl.clf()
        pl.plot(declination,shadowing)
        pl.plot(declination,shadowing,'o')
        pl.xlabel("Declination (Degree)")
        pl.ylabel("Shadowing (%)")
        pl.savefig("Shadowing-Declination.png")
        pl.show()
    

class figures:
    """
    Class to produce different plots
    """
    
    def __init__(self):
        
        pass
    
    
    def arrayConfig(self,fileCfg,xmin,xmax,ymin,ymax,title="",auto=False,savefig = None):
        """
        fileCfg: configuration file for the Array
        xmin,xmax,ymin,ymax:  plot ranges
        title: title in the plot
        """
    
        xx=[]
        yy=[]
        zz=[]
        name=[]
        
        f=open(fileCfg)
        
        dump=f.readline()
        res=dump.split()
    
    
        while dump !="":
            if res[0] != "#":
                xx.append(res[0])
                yy.append(res[1])
                zz.append(res[2])
                name.append(res[4])
            dump=f.readline() 
            res=dump.split()  
            print res
            
        f.close()
        
        xa=np.array(xx,np.float)
        ya=np.array(yy,np.float)
        
        
        xmean=xa.mean()
        ymean=ya.mean()
        
        xa-=xmean
        ya-=ymean
        
        if auto :
            xmin = 1.1 * min(xa)
            xmax = 1.1 * max(xa)
            ymin = 1.1 * min(ya)
            ymax = 1.1 * max(ya)
        
        pl.clf()
        pl.plot(xa,ya,"ro")
        pl.xlabel("X (meter)")
        pl.ylabel("Y (meter)")
        pl.xlim(xmin,xmax)
        pl.ylim(ymin,ymax)
               
        index=0
        for ant in name:
            pl.text(xa[index],ya[index],ant)
            index+=1
        pl.text(xmin*0.9,ymax*0.9,title)
        pl.show()


class report:
    """
    Meta-class to create all the plots and reports.
    """
    
    def __init__(self,directoryFile,reportType):
        
        self.directoryFile=directoryFile
        self.reportType-reportType
        self.band={'3':'100GHz','6':'230GHz','7':'345GHz','9':'675GHz'}
        
        
    def  create(self,filename):
        """
        Create the report,text+plots+simuls, according to the reportType
        filename : Antenna Configuration file
        """
        
        
        pass
   
        
    
    
########################Main program####################################
if __name__=="__main__":
    " main program"      
        
    a=Shadowing()
    dec,shad=a.curveShadowingDeclination('C32-cpt.cfg','2h',-75,25,2)
    a.plotShadowingDeclination(dec,shad)
    out=open("shadowingDec-2h.pickle","wb")
    data={'dec':dec,'shad':shad}
    pickle.dump(data,out)
    out.close()
 
 

    #f=figures()   
    #f.arrayConfig("Pads.cfg",-1000,1000,-1000,1000)
    #f.arrayConfig("extendedCycle0.cfg",-210,270,-240,240)
    
    

    






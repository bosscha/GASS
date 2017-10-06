#!/usr/bin/python

"""
Class/methods to analyze the properties of a dirty beam (sidelobes, size)



HISTORY:
    2011.03.12:
        - first shot
        - read the image of a beam (.psf)
        - compute the radial distribution of the maximum (absolute) of an image (typically .psf) 

    2011.04.16:
        - plot the maxima of the radial level of a beam
        
    20013.01.23:
        - updating beam to work with CASA 4.0.0

    2014.10.24:
        - change the name to beamAnalysis
        - add a class to analyze properties of the beam for an array configuration
        


RUN:
In casapy, you can execute:
> execfile('beamAnalysis')

or 
> from beamAnalysis import *

Example of use:

> b=beam("3c288-60-snap.psf")
> b.plotRadialMaxLevel(0.,10.,100)  

"""

__authors__="ALMA: SL"
__version__="0.3.0@2014.10.24"


import numpy as  np
from  math import * 
from casa import image as ia
import pylab as pl
from  casa import componentlist as cl
from casa import simobserve, simanalyze
import os
import inspect
import casa


RAD2ARCSEC=3600.*180./pi

class beam:
    
    
    def __init__(self,beamImage):
        
        self.beamImage=beamImage
        
        self.__setImage()
        
        
    def __setImage(self):
        "Set the parameters of the image"
        
        ia.open(self.beamImage)
        hdr = ia.summary()
        
        print hdr
        
        ## Image Size
        self.Nx=hdr['shape'][0]
        self.Ny=hdr['shape'][1]
        
        ## Pixel reference
        self.refx=hdr['refpix'][0]
        self.refy=hdr['refpix'][0]
        
        ##Increment (radians)
        self.dx=hdr['incr'][0]
        self.dy=hdr['incr'][1]
        self.dx*=RAD2ARCSEC
        self.dy*=RAD2ARCSEC
        
        ia.close()
        

    def radialMaxLevel(self,rmin,rmax,Nbin):
        "Compute the maximum (absolute) level for the image and return the array of radius (arcsec) and values"
        
        dr=(rmax-rmin)/Nbin
        radius=np.arange(Nbin)*dr+rmin
        radMaxLevel=np.zeros(Nbin)
       
        
        
        ## loop over the pixels
        
        ia.open(self.beamImage)
        
        for i in range(self.Nx):
            for j in range(self.Ny):
                xx=(i-self.refx)*self.dx
                yy=(j-self.refy)*self.dy
                rr=sqrt(xx*xx+yy*yy)
                indexR=floor((rr-rmin)/dr)
                
                if indexR >=0 and indexR < Nbin:
                    val=ia.pixelvalue([i,j])['value']['value']
                    if abs(val) > abs(radMaxLevel[indexR]):
                        radMaxLevel[indexR]=val
        
        ia.close()
        
        return(radius,radMaxLevel)
    
    
    def plotRadialMaxLevel(self,rmin,rmax,Nbin):
        """
        Plot the maximum radial levels of a beam and  a 5 and 10 % line
        rmin,rmax: min and max angular size (arcsec)
        """
        
        rr,val=self.radialMaxLevel(rmin,rmax,Nbin)
        
        pl.plot(rr,100*val)
        
        ymin=min(100*val)
        
        pl.xlim([rmin,rmax])
        pl.ylim([ymin*1.05,100.])
        
        pl.plot([rmin,rmax],[5.,5.],'r--')
        pl.plot([rmin,rmax],[0.,0.],'k-')
        pl.plot([rmin,rmax],[-5.,-5.],'r--')
        pl.text(rmax*0.85,6.,"5% level")
        
        pl.xlabel('radius (arcsec)')
        pl.ylabel("Beam max. level (%)")
                
        pl.show()
                
                
class arrayConfigurationBeam:
    
    def __init__(self, fileCfg):
        "CASA antenna configuration file"
        
        self.fileCfg = fileCfg
        
        
    def curveBeamResolutiongDeclination(self, trackDuration, declinationMin,declinationMax,decStep):
        """
        Compute the minor,Major axis vs. Declination for a given configuration
            trackDuration: time duration, e.g. "6h"
            declinationMin,declinationMax: range of declination
            decStep: step for the declination, needs to be an integer...
            
        Output : arrays [declination] [minor] [major]
        """
              
        antennaCfg = self.fileCfg
        print antennaCfg
        
        projectName = "tempCurveBeamDeclination"
        nStep = floor((declinationMax-declinationMin)/decStep)
        decArr = np.arange(nStep)*decStep+declinationMin
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
            
                    
            sim = inspect.getargspec(casa.simobserve)
            
            ind = 0
            for arg in sim.args:
                casa.simobserve.__setattr__(arg, sim.defaults[ind])
                ind += 1
                                      
            casa.simobserve.__setattr__('project',projectName+"%d"%(index))
            casa.simobserve.__setattr__('complist',projectName+"%d.cl"%(index))
            casa.simobserve.__setattr__('compwidth','8GHz')
            casa.simobserve.__setattr__('antennalist',antennaCfg)
            casa.simobserve.__setattr__('totaltime',trackDuration)
            casa.simobserve.__setattr__('integration','10s')
            
            casa.simobserve(antennalist = 'best-O-1_step2_88.cfg' , complist = 'tempCurveBeamDeclination0.cl')
                    
                       
                       
            simanalyze(
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
            imageName = projectName+"%d/"%(index)+projectName+"%d.%s.image"%(index,antennaCfg.split('.')[0])
            
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
            las100[index] = RAD2ARCSEC * 0.0017987547479999 / duvMin      ## 0.6 * Lambda / BL_min
            
            print "mininum BL:",duvMin
            print "LAS (100 GHz):",las100[index]
            
            
            
            index+=1
            
        
        # fout=open("curveBeamDec.last","wb")
        # data={'Dec':decArr,'Shadow':fractionShadowing}
        # pickle.dump(data,fout)  
        # fout.close()  
        
        return(decArr,minorAxis,majorAxis,las100)
        
   
###################Main progra*m########################
if __name__=="__main__":
    " main program"
    
    print "BeamAnalysis.py external task (still in testing...)"
    
    

    
    
  
        
        
#!/usr/bin/python

"""
Class to analyze UVW visibilities  within CASA

HISTORY:
    2011.03.12:
        - first shot
        - class UVW to analyze the uv-coverage in a MS
        
    2011.03.13:
        - radial density of UV distance.
        - plot of the radial density (normalized to the number of visibilities not 1!!)
        - plot of the UV coverage
        
    2011.03.31:
        - find the shadowed antennas in the MS (from Remy Indebetouw)
        
    2011.04.04:
        - Azimuth Density: computation+plot
        
    2011.04.16:
        - Add noShadowing method to select good antennas (still bugged!!)
        
    2011.05.10:
        - minor changes in the  plot functions
        
    2011.05.15:
        - add Amplitude in the UVW class
        - Add plot of a list ms for the Amp. vs. Ruv. Still bugged...
        
    2011.05.17:
        - Modify the plot UVDist_Amp
        
    2011.11.07:
        - modify the import to be used in casa properly. No need to do a execfile in principle.
        
    2011.11.08:
        - fix shadowing method
        
    2011.11.10:
        - add figfile in plotUVcoverage
        
    2012.03.08:
        - add a class  setUVW to manipulate the visibilities in a MS
        
    2013.02.01:
        - quartile (min,25%,50%,75%,max) of the uv-distance distribution
        
    2013.08.07:
        - add Norm option in radial density  to count Nbl only.
        
    2017.07.21:
        - adapting to CASA version >= 5.0.0
        
        
RUN:

In casapy, to start execute:
> execfile('UVW.py')

Examples of use:

## Radial Density Plot
> a=UVW("3c288.ms")
> a.plotRadialDensity(0.,130.,100)   

## Extrema (meter) of the UV distance
> a=UVW("3c288-30.ms")
> ruv=a.distUV()
> print min(ruv),max(ruv)

## UV coverage plot
>a=UVW("3c288-30.ms")
>a.plotUVCoverage(-130.,130.,-130.,130.,"V. vs. U (Dec=-30)")

#### To add shadowing (red) on the plot
>a.plotUVCoverage(-130.,130.,-130.,130.,"V. vs. U (Dec=-30)",shadow=True,saveFig=True,filefig='toto.png')


use:
sys.path.insert(0,'/home/stephane/git/ALMA/ALMA/ArrayConfiguration/')
"""

__authors__="ALMA : SL"
__version__="0.4.0@2017.07.21"




import scipy as sci
import numpy as np
import math

from casa import *
from taskinit import *

# from casa import table as tb
import pylab as  pl

RAD2DEGREE = 180. / math.pi
DEGREE2RAD = math.pi / 180.

# ms = casac.ms()

class UVW:
    
    
    def __init__(self,vis):
        "vis is a MS file"
        
        self.vis=vis
        self.__readUVW()
        
        
    def __readUVW(self):
        """Read the UVW data from the MS. Works fine for a simulation (simdata). 
        To Be Checked for a real obs. with different targets."""
        
        ms.open(self.vis)
        uvw = ms.getdata(['amplitude','u','v','w'])
        ms.done()
        
        tb.open(self.vis)
        self.f=tb.getcol('FLAG_ROW')
        self.ff=tb.getcol('FLAG')
        tb.done()
        
        self.u=uvw['u']
        self.v=uvw['v']
        self.w=uvw['w']
        
        self.amp=uvw['amplitude']
           
            

        
    def distUV(self,noShadow= False):
        "Return the module of the UV vector"
        
        if not noShadow:
            rr2=self.u*self.u+self.v*self.v
            duv=np.zeros(len(rr2))
        else:
            nSh = self.noShadowing()
            rr2=self.u[nSh]*self.u[nSh]+self.v[nSh]*self.v[nSh]
            duv=np.zeros(len(rr2))
            
        for i in range(len(rr2)):
            duv[i]=math.sqrt(rr2[i])
        
        return(duv)
    
    
    
    def quartileUV(self):
        "Return the quartile (min,25%,50%,75%,max) of the uv-distance in the MS"
        
        duv = self.distUV()
        nData = len(duv)
        
        duvSort = np.sort(duv)
        
        minBL = duvSort[0]
        maxBL = duvSort[nData-1]
        Q1 = duvSort[int(0.25*(nData-1))]
        Q2 = duvSort[int(0.50*(nData-1))]
        Q3 = duvSort[int(0.75*(nData-1))]
        
        return([minBL,Q1,Q2,Q3,maxBL])
        
    
    def radialDensity(self,minRadius,maxRadius,nBin, norm = True):
        "Return the radius (meter) and the radial density array"
        
        duv=self.distUV()
        nData=len(duv)
        
        dr=(maxRadius-minRadius)/nBin
        
        radius=np.zeros(nBin)
        rhoUV=np.zeros(nBin)
        
        #Building the density
        for value in duv:
            indexFloating=(value-minRadius)/dr
            index=int(indexFloating+0.5)
            
            if index >=0 and index < nBin:
                rhoUV[index] +=1.
                
        #Normalization 
        for i in range(nBin):
            r=minRadius+(i+0.5)*dr
            dS=2.*math.pi*r*dr
            radius[i]=r
            if norm :
                rhoUV[i]=rhoUV[i] / dS
            else :
                ds = 1.
                
                       
        return(radius,rhoUV)

    def azimuthDensity(self,nBin):
        "Compute the Azimuth Density between 0 and 180 degrees"
        
        nData=len(self.u)
        dangle=180./nBin


        
        angle=np.arange(nBin)*dangle
        rhoAzimuth=np.zeros(nBin)
        
        
        for i in range(nData):
            
            angleUV=math.atan2(self.v[i],self.u[i])*RAD2DEGREE
            indexAzimuthFloating=angleUV/dangle
            index=math.floor(indexAzimuthFloating+0.5)
  
            if index>=0 and index < nBin:
                rhoAzimuth[index]+=1.
  
            angleUV=math.atan2(-self.v[i],-self.u[i])*RAD2DEGREE
            indexAzimuthFloating=angleUV/dangle
            index=math.floor(indexAzimuthFloating+0.5)
  
            if index>=0 and index < nBin:
                rhoAzimuth[index]+=1.
  
  
            
        return(angle,rhoAzimuth)
            
            
        
        

    def shadowing(self):
        "Select the shadowed antennas from the FLAG column and return the index of the shadowed measurement and the percentage of shadowing "
        
        indexFlag=pl.concatenate((pl.where(self.f==1)[0],pl.where(self.ff[0,0,])[0],pl.where(self.ff[1,0,])[0]))
        indexNoFlag=pl.concatenate((pl.where(self.f==0)[0],pl.where(self.ff[0,0,]==False)[0],pl.where(self.ff[1,0,]==False)[0]))
        
        Ntot=len(indexFlag)+len(indexNoFlag)
        fractionShadow=100.*len(indexFlag)/Ntot
        
        return(indexFlag,fractionShadow)
        
 
    def noShadowing(self):
        "Select the non flagged (non shadowed) antennas from the FLAG and the index... Still bugged !!"
        
        indexNoFlag=pl.concatenate((pl.where(self.f==0)[0],pl.where(self.ff[0,0,]==False)[0],pl.where(self.ff[1,0,]==False)[0]))
        
        return(indexNoFlag)
 
            
            
 

    def plotRadialDensity(self,minRadius,maxRadius,nBin,savefig=False,figfile='toto.png'):
        "Plot the UV radial density (Normalized to the number of visibilities)"
        
        rr, rho = self.radialDensity(minRadius,maxRadius,nBin)
        
        pl.xlabel(r'$R_{UV} (meter)$',fontsize=15)
        pl.ylabel(r'$\rho (R)$',fontsize=15)
        # pl.xlim(0.,140.)
        # pl.ylim(0.,10.)
        #pl.grid(True)
        pl.plot(rr,rho,"b-")
        
        if savefig:
            pl.savefig(figfile)
        else: 
            pl.show()
        pl.clf()
        
        
        
    def plotAzimuthDensity(self,nBin,title,savefig=False,figfile='toto.png'):
        "Plot the Azimuth Density normalized to the 0 degree value"
        
        ang,dens=self.azimuthDensity(nBin)
        dens/=dens[0]
        
        ymax=max(dens)
        
        pl.clf()
        pl.plot(ang,dens)
        pl.xlim([0.,180.])
        pl.ylim([0.,ymax*1.05])
        
        pl.xlabel("azimuth (Degree)",fontsize=15)
        pl.ylabel(r'$\rho (\theta)$',fontsize=15)
        pl.title(title)
        
        if savefig:
            pl.savefig(figfile)
        else: 
            pl.show()
        pl.clf()
        
        
        
        
    def plotUVCoverage(self,umin,umax,vmin,vmax,title,shadow=False,saveFig=False,figfile='toto.png'):
        "UV coverage plot with the UV extrema. If shadow=True, will plot the shadowing in red"
        
        l=pl.plot(self.u,self.v,"k.")
        pl.setp(l, 'markersize', 0.4)
        l=pl.plot(-self.u,-self.v,"k.")
        pl.setp(l, 'markersize', 0.4)
        
        if shadow:
            indexShadow,fracShadow=self.shadowing()
            l=pl.plot(self.u[indexShadow],self.v[indexShadow],"r.")
            pl.setp(l, 'markersize', 0.4)
            l=pl.plot(-self.u[indexShadow],-self.v[indexShadow],"r.")
            pl.setp(l, 'markersize', 0.4)
            
            print("Fraction of shadowing: %4.2f %%"%(fracShadow))
        
                   
        pl.xlim(umin,umax)
        pl.ylim(vmin,vmax)
        pl.title(title,fontsize=15)
        pl.xlabel("U (meter)",fontsize=15)
        pl.ylabel("V (meter)",fontsize=15)
        pl.xticks(fontsize=15)
        pl.yticks(fontsize=15)
        
        if saveFig:
            pl.savefig(figfile)
        else :
            pl.show()
            
        pl.clf()

def plotUVDist_Amp(listOfMs1,listOfMs2,label1,label2,title):
    """
    plot the Amplitude vs. the UV distance for a list of .ms
    """
    
        
    pl.xlabel(r'$R_{UV} (meter)$',fontsize=15)
    pl.ylabel("Amplitude",fontsize=15)
    pl.text(title[0],title[1],title[2])
    
    index=0    
    for msName in listOfMs1:
        
        print msName
        UVWms=UVW(msName)
        dUV=UVWms.distUV()
        
        print dUV
        print UVWms.amp[0][0]
        pl.plot(dUV,UVWms.amp[0][0],"k.",markersize=0.5)
        pl.text(label1[index][0],label1[index][1],label1[index][2])
        index+=1
  
        
    index=0
    for msName in listOfMs2:
        
        print msName
        UVWms=UVW(msName)
        dUV=UVWms.distUV()
        
        print dUV
        print UVWms.amp[0][0]
        pl.plot(dUV,UVWms.amp[0][0],"r.",markersize=0.5)
        pl.text(label2[index][0],label2[index][1],label2[index][2],color='red')
        index+=1
        
    pl.show()
    
    
class setUVW():
    "To manipulate the visibilities of a MS"
    
    def __init__(self,visName):
        
        self.visName = visName
        
        
    def addPhasePerAntenna(self,phase,antennaId,dataType="DATA"):
        "Add a phase (degree) per antennaId. Careful with the order !"
        
        tb.open(self.visName,nomodify=False)
        
        ant1 = tb.getcol("ANTENNA1")
        ant2 = tb.getcol("ANTENNA2")
        
        data = tb.getcol(dataType)

        phase = DEGREE2RAD*phase

        print "Phase rotation: %f"%(phase)

        nDim = data.shape
        
        for i in range(nDim[0]):
            for j in range(nDim[1]):
                for k in range(nDim[2]):
                    antenna1 = ant1[k]
                    antenna2 = ant2[k]
                    
                    
                    if antenna1 == antennaId:
                        value = data[i][j][k]
                        real = value.real
                        image = value.imag
                        cosRot = math.cos(phase)
                        sinRot = math.sin(phase)
                        
                        xRot = cosRot*real-sinRot*image
                        yRot = sinRot*real+sinRot*image
                        
                        data[i][j][k] = complex(xRot,yRot)
                        
                        # print "Amplitude:%f"%(math.sqrt(xRot*xRot+yRot*yRot))
 
                        
                    if antenna2 == antennaId:
                        value = data[i][j][k]
                        real = value.real
                        image = value.imag
                        cosRot = math.cos(-phase)
                        sinRot = math.sin(-phase)
                        
                        xRot = cosRot*real-sinRot*image
                        yRot = sinRot*real+sinRot*image
                        
                        data[i][j][k] = complex(xRot,yRot)
                        # print "Amplitude:%f"%(math.sqrt(xRot*xRot+yRot*yRot))
                        
                    
        tb.putcol(dataType,data)
        tb.close()
        
        
            
            
    
    
########################Main program####################################
if __name__=="__main__":
    " main program"
    
    print "UVW.py external task (still in testing...)"
    
    #a=UVW("1mm.d-20.9asec.ms")
    #a.plotAzimuthDensity(50,"test")
    
    
    ####### 6h integration in the extended configuration
    #list1=['100GHz6h.1asec.ms','100GHz6h.4asec.ms','100GHz6h.7asec.ms']
    #list2=['230GHz6h.1asec.ms','230GHz6h.4asec.ms','230GHz6h.7asec.ms']
    #label1=[[173,0.78,'1"'],[88,0.36,'4"'],[65,0.2,'7"']]
    #label2=[[75,0.65,'1"'],[20,0.3,'4"'],[18,0.05,'7"']]   
    #plotUVDist_Amp(list1,list2,label1,label2,title=[340,0.91,"6h integration"])
    
    
    ####### 30m integration in the extended configuration
    #list1=['100GHz30m.1asec.ms','100GHz30m.4asec.ms','100GHz30m.7asec.ms']
    #list2=['230GHz30m.1asec.ms','230GHz30m.4asec.ms','230GHz30m.7asec.ms']
    #label1=[[173,0.78,'1"'],[88,0.36,'4"'],[65,0.2,'7"']]
    #label2=[[75,0.65,'1"'],[21,0.3,'4"'],[26,0.022,'7"']]   
    #plotUVDist_Amp(list1,list2,label1,label2,title=[340,0.91,"30m integration"])
        
        
        
    
    
    
        
        
        
    
    
    
            
        
        
        
        
        
        
        
        
        
        


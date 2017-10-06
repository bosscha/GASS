#!/usr/bin/python


"""
Script to return the Cycle 2 arrays for a given LAS, Angular Resolution


HISTORY:
    2012.10.11:
        - first shot
        - Cycle 1 setup
        
    2012.11.20:
        - Adding resolution, LAS for different declination
        
    2012.12.27:
        - adding a filter if no configuration is found to multiply by a fudge factor
        
        
    2013.03.04:
        - Adding geometric average for the spatial resolution  and test on twice the spatial resolution
        
    2013.03.05:
        - fixing bugs
        - adding extra information (THB, array)
        
        
    2013.03.11:
        - Removing the condition about AR/2. for the array.
        
        
    2013.03.18:
        - Adding point source option
        - putting the AR/2. if it is not a point source
        
        
    2013.03.20:
        - changing the PS option to LAS = 0
        
    2013.05.02:
        - changing slightly the conditions of acceptance (>= instead of >)
        
    2013.05.03:
        - print version
        - try increasing fudge factor until it gets a solution

    2013.05.10:
        - Add silentRun for P2G (FG)
        
    2013.12.13:
        - Name change for Cycle2 and update of the pickle.
        - Update of the finder (multi-configuration)
        
    2013.12.16:
        - new find_array3 for the multi-configuration
    
    
    2014.01.28:
        - fix the fudge relaxation for the resolution.
        
        
    2014.05.22:
        - New algorithm to deal with minAr, maxAr
        
RUN:
Input RES (arcsec) LAS (arcsec) FREQ (GHz) Declination (Degree) Y/N (ACA)
> python arrayResolution2.py 0.2 2.5 640. -53 Y PS

"""

__author__ = "ALMA : SL, AL, FG"
__version__ = "0.6.0@2014.05.22"


import sys, pickle
import math


### ALMA
LATITUDE_ALMA = -23.03
DEG2RAD = math.pi/180.

class arrayRes:
    
    def __init__(self, arguments):
        
        self.LAS = [26.1,26.3,18.0,18.0,14.4,9.1,9.1]           
        self.LASA = [44.0,44.0,44.0,44.0,14.4,9.1,9.1]          
        self.LAST = [390.0,390.0,390.0,390.0,14.4,9.,9.1]      
        self.res = [3.73,2.04,1.40,1.11,0.75,0.57,0.41]              
            
        self.frequencyReference = 100.
        self.lasObs = 0.        
        self. resObs = 0.
        self.freqObs = 0.
        self.pointSource = False
        self.silent = True
        
        self.args = arguments 
                      
        self.array = {0:"C34-1",1:"C34-2",2:"C34-3",3:"C34-4",4:"C34-5",5:"C34-6",6:"C34-7"}
        self.read_cycle2()
        
        

    def set_las(self,las):
        "Set the LAS of the observation"
        
        self.lasObs = las   
        
    def set_res(self,res):
        "Set the angular resolution of the observation"
        
        self.resObs = res 
        
    def set_frequency(self,freq):
        "Set the frequency of the observation"
        
        # if ((freq>=64.)&(freq<=116.)): freq = 100.
        # if ((freq>=211.)&(freq<=275.)): freq = 230.
        # if ((freq>=275.)&(freq<=373.)): freq = 345.
        # if ((freq>=602.)&(freq<=720.)): freq = 675.
        
        self.freqObs = freq
                
    
    def set_declination(self,declination):
        "Set the representative declination of the observation"
        
        self.declination = declination
        
        
    def set_aca(self,aca):
        "Set the frequency of the observation"
        
        self.acaUse = aca

    
    def set_pointSource(self, isPS):
        "Set True if point source"
        
        self.pointSource = isPS
    
    def read_cycle2(self,directory = './'):
        "Read the detailed Resolution, LAS for the Cycle 1 configuration"
        
        
        self.data = []
        
        f = open(directory+'Resolution-C34-1.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C34-2.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C34-3.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C34-4.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C34-5.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C34-6.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C34-7.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        
        ### ACA ####
        f = open(directory+'Resolution-ACA-std.pickle')
        self.aca = pickle.load(f)
        f.close()
        
        
        

    def find_array(self):
        "Find the array with the obs. input"
        
        TP='N'
        arrayMatch = []
        scalingFrequency = self.frequencyReference / self.freqObs

        if (self.acaUse == 'Y'):
            self.LAS=self.LAST
            if self.lasObs / scalingFrequency > self.LASA[1]:
                TP='Y'
                
        for arr in self.array :
            
            if self.silent:
                print  self.LAS[arr] * scalingFrequency, self.res[arr] * scalingFrequency
            
            if self.LAS[arr] * scalingFrequency >= self.lasObs and self.res[arr] * scalingFrequency <= self.resObs:
                arrayMatch.append(self.array[arr])
            else:
                arrayMatch.append("")
        return arrayMatch,TP
            
            
            
    def find_array2(self,verbose = False):
        "Find the array with the obs. input using the precise resolution, LAS..."

        TP = 'N'
        scalingFrequency = self.frequencyReference / self.freqObs


        nData = len(self.data[0][0])
        decMin = self.data[0][0][0]
        decMax = self.data[0][0][nData-1]
        deltaDec = (decMax-decMin)/nData
        
        index = int(math.floor(((self.declination-decMin) / deltaDec)))

        # print index
        
        ### No ACA
        arrayMatch = []
        
        for arr in self.array :   

            lasArr = self.data[arr][3][index]
            resArr = math.sqrt(self.data[arr][1][index] * self.data[arr][2][index])
            
            lasFreqArr = lasArr * scalingFrequency
            spatialResolutionArr = resArr * scalingFrequency
            
            res_thb = self.res[arr]*scalingFrequency
            las_thb = self.LAS[arr]*scalingFrequency
            elevation_factor = abs(1./math.sin(DEG2RAD*(90.-LATITUDE_ALMA+self.declination)))
            
            res_estimated = math.sqrt(res_thb*res_thb*elevation_factor)
            las_estimated = math.sqrt(las_thb*las_thb*elevation_factor)

            if self.silent:
                if(verbose):
                    print("# Array: %s, LAS: %5.2f, RES: %5.2f"%(self.array[arr],lasFreqArr, spatialResolutionArr )) 
                    print("#          THB: LAS: %5.2f, RES: %5.2f")%(las_estimated,res_estimated)
                    # print("#")  
                 
            if self.pointSource:
                if lasFreqArr >= self.lasObs and self.resObs >= spatialResolutionArr  :
                    arrayMatch.append(self.array[arr])
                else:
                    arrayMatch.append("")
                    
            else :
                if lasFreqArr >= self.lasObs and self.resObs >= spatialResolutionArr and spatialResolutionArr  >= self.resObs / 2.  :
                    arrayMatch.append(self.array[arr])
                else:
                    arrayMatch.append("")

        ### ACA  used 
        if (self.acaUse == 'Y'):
            arrayMatch = []
            
            for arr in self.array:
                    resArr = math.sqrt(self.data[arr][1][index] * self.data[arr][2][index])
                    spatialResolutionArr = resArr*scalingFrequency
                    
                    ## 
                    
                    if self.pointSource:
                        if self.resObs > spatialResolutionArr   and arr < 4:
                            arrayMatch.append(self.array[arr])
                        else:
                            arrayMatch.append("")
                            
                    else :
                        if self.resObs >= spatialResolutionArr   and spatialResolutionArr  >= self.resObs / 2. and arr < 4:
                            arrayMatch.append(self.array[arr])
                        else:
                            arrayMatch.append("") 
                            
            
            lasACA = self.aca[3][index]
            if lasACA*scalingFrequency <= self.lasObs:
                TP = 'Y'
            
            
        return arrayMatch, TP
      
      
    def find_array3(self,verbose = False):
        "Find the array with the obs. input using the precise resolution, LAS.... It takes into account a multi-configuration"

        TP = 'N'
        scalingFrequency = self.frequencyReference / self.freqObs
	

        nData = len(self.data[0][0])
        decMin = self.data[0][0][0]
        decMax = self.data[0][0][nData-1]
        deltaDec = (decMax-decMin)/nData
        
        index = int(math.floor(((self.declination-decMin) / deltaDec)))

        # Cycle 2 Match Array
        
        matchArrayCycle2 = {3:0,4:1,5:2,6:2}
        
        ###
        arrayMatchRes = []
        arrayMatchLAS = []
        
        lasFreqArrAll = []
        resFreqArrAll = []
       
       
       
	for arr in self.array :   
	    arrayMatchRes.append("")
	    arrayMatchLAS.append("")
	    
            lasArr = self.data[arr][3][index]
            resArr = math.sqrt(self.data[arr][1][index] * self.data[arr][2][index])
            
            lasFreqArr = lasArr * scalingFrequency
            spatialResolutionArr = resArr * scalingFrequency
            
            lasFreqArrAll.append(lasFreqArr)
            resFreqArrAll.append(spatialResolutionArr)
            
            res_thb = self.res[arr]*scalingFrequency
            las_thb = self.LAS[arr]*scalingFrequency
            elevation_factor = abs(1./ math.sin(DEG2RAD*(90.-LATITUDE_ALMA+self.declination)))
            
            res_estimated = math.sqrt(res_thb*res_thb*elevation_factor)
            las_estimated = math.sqrt(las_thb*las_thb*elevation_factor)

            if self.silent:
                if(verbose):
                    print("# Array: %s, LAS: %5.2f, RES: %5.2f"%(self.array[arr],lasFreqArr, spatialResolutionArr )) 
                    print("#          THB: LAS: %5.2f, RES: %5.2f")%(las_estimated,res_estimated)
                    # print("#")  
                          
	########################### Comparison #######################
         
  	notFound = True   
  	notFoundLAS = True
            
        for arr in self.array :   
            
            lasFreqArr = lasFreqArrAll[arr]
            spatialResolutionArr = resFreqArrAll[arr]
            
 
	    if self.pointSource:
                if self.resObs >= spatialResolutionArr  :
                    arrayMatchRes[arr] = self.array[arr]
                    notFound = False
                                 
            else :
	      
		
                if  self.resObs >= spatialResolutionArr and spatialResolutionArr  >= self.resObs / 2.  :
                    arrayMatchRes[arr] = self.array[arr]
                    notFound = False
                    
                    if lasFreqArr <= self.lasObs and arr > 2:
		      arrayMatchLAS[matchArrayCycle2[arr]] = self.array[matchArrayCycle2[arr]]
		      
		      if lasFreqArrAll[matchArrayCycle2[arr]] <= self.lasObs and matchArrayCycle2[arr] > 0:
			for i in range(0,matchArrayCycle2[arr]):
			  if lasFreqArrAll[i] >= self.lasObs :
			    arrayMatchLAS[i] = self.array[i]
			    notFoundLAS = False
			   
		   
                    
                
        ### ACA  used  ###############
        
        if (self.acaUse == 'Y'):
	  
	  
	  arrayMatchRes = []
	  arrayMatchLAS = []
       
	  for arr in self.array :   
	    arrayMatchRes.append("")
	    arrayMatchLAS.append("")
           
            
          for arr in self.array:
            spatialResolutionArr = resFreqArrAll[arr]                                            
               
            if  self.resObs >= spatialResolutionArr and spatialResolutionArr  >= self.resObs / 2.  :
	      arrayMatchRes[arr] = self.array[arr]
              notFound = False
              
              if  arr > 2:
		arrayMatchLAS[matchArrayCycle2[arr]] = self.array[matchArrayCycle2[arr]]
 
 
            lasACA = self.aca[3][index]
            if lasACA*scalingFrequency <= self.lasObs:
                TP = 'Y'
            
            
        return [arrayMatchRes,arrayMatchLAS] , TP , notFound, notFoundLAS
    
    
    
    
    def find_putAR(self, verbose = False):
        " Return (LAS) minAr, maxAR per configuration (not the configuration)"
        
        pass 
      
  ########################################################################3       
    def silentRun(self):
        self.silent = False
        
    def run(self):
        "Run the matching array"
        
        TP="N"

        self.set_res(float(self.args[1]))
        self.set_las(float(self.args[2]))
        self.set_frequency(float(self.args[3]))
        self.set_declination(float(self.args[4]))
        self.set_aca((self.args[5]))
        
        if self.lasObs == 0.:
                self.set_pointSource(True)
        
        
        strOut = "### arrayResolution2 \n"
        strOut += "### Version: %s \n"%(__version__)
        strOut += "### Input \n"
        if self.pointSource:
            strOut += "# Point Source ! \n"
        strOut += "# Spatial Resolution: %s \n"%(self.args[1])
        strOut += "# LAS: %s \n"%(self.args[2])
        strOut += "# Frequency: %s GHz \n"%(self.args[3])
        strOut += "# Declination: %s \n"%(self.args[4])
        strOut += "# ACA (Y/N): %s \n\n"%(self.args[5])
        strOut += "### Output (target frequency) \n"
        strOut += "### Using CASA simulation with natural weighting (slightly different of THB)"

        if self.silent:        
            print(strOut)
                                                

        notFound = True
        maxTry = 100
        nTry = 1
        deltaFudgeFactor = 0.05
        fudgeFactor = 1.0
        
        res , TP , notFound , notFoundLAS = self.find_array3(verbose = True)
           
        
            
        while (notFound and nTry < maxTry and notFoundLAS and self.acaUse == 'N' ):
            
            nTry += 1
            notFound = False
            fudgeFactor += deltaFudgeFactor
            self.resObs *= fudgeFactor
            self.lasObs /= fudgeFactor
            res , TP , notFound , notFoundLAS = self.find_array3()
            
        
        while (notFound and nTry < maxTry ):
            
            nTry += 1
            notFound = False
            fudgeFactor += deltaFudgeFactor
            self.resObs *= fudgeFactor
            res , TP , notFound , notFoundLAS = self.find_array3()
            
        
        
            
        if nTry > 1 :
            if self.silent:
                print("# No array configuration found, fudge factor applied (Tol = %3.0f %%)"%((fudgeFactor-1.)*100.))
        
        if notFound and nTry > 1:
            if self.silent:
                print("# No array configuration found, even with fudge factor, problem ...")
            

        strOutRes = ""
        strOutLAS = ""
        pcomR = ""
        pcomL = ""

        if self.silent:        
            print ""
            print("### Results  - AR - LAS")

        if notFound :
            if self.silent:
                print ",,,,,"
        else:
            for s in res[0]:
                strOutRes += pcomR+s
                pcomR = ","
	    strOutRes += ","
	    strOutRes += TP
	    
	    for s in res[1]:
                strOutLAS += pcomL+s
                pcomL = ","
	    strOutLAS += ","
	    strOutLAS += TP

        if self.silent:        
            print strOutLAS
            print strOutRes
        
        return strOutRes , strOutLAS
        
#====== Standalone program =========================
if __name__=="__main__":
       
    
    arg = sys.argv

    if len(arg) < 6:
        print "Arguments missing \n"
        print "The correct syntax is:"
        print "python arrayResolution.py RES (arcsec) LAS (arcsec) FREQ (GHz) Declination (Degree) Y/N (ACA) \n"
        print "Example:"
        print "python arrayResolution2.py 0.2 0. 640. -53 Y    ## if LAS = 0. assumes a point source"
       
        
    else :
        a = arrayRes(arg)
        a.run()
        
        
        
    

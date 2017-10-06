#!/usr/bin/python


"""
Script to return the Cycle 1 arrays for a given LAS, Angular Resolution


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
        
RUN:
Input RES (arcsec) LAS (arcsec) FREQ (GHz) Declination (Degree) Y/N (ACA)
> python arrayResolution.py 0.2 2.5 640. -53 Y PS

"""

__author__ = "ALMA : SL, AL, FG"
__version__ = "0.4.0@2013.05.20"


import sys, pickle
import math


### ALMA
LATITUDE_ALMA = -23.03
DEG2RAD = math.pi/180.

class arrayRes:
    
    def __init__(self, arguments):
        
        self.LAS = [25.0,25.0,17.0,17.0,14.0,8.6]
        self.LASA = [44.0,44.0,44.0,44.0,14.0,8.6]
        self.LAST = [390.0,390.0,390.0,390.0,14.0,8.6]
        self.res = [3.7,2.0,1.4,1.1,0.75,0.57]
            
        self.frequencyReference = 100.
        self.lasObs = 0.        
        self. resObs = 0.
        self.freqObs = 0.
        self.pointSource = False
        self.silent = True
        
        self.args = arguments 
                      
        self.array = {0:"C32-1",1:"C32-2",2:"C32-3",3:"C32-4",4:"C32-5",5:"C32-6"}
        self.read_cycle1()
        
        

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
    
    def read_cycle1(self,directory = './'):
        "Read the detailed Resolution, LAS for the Cycle 1 configuration"
        
        
        self.data = []
        
        f = open(directory+'Resolution-C32-1.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C32-2.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C32-3.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C32-4.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C32-5.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        f = open(directory+'Resolution-C32-6.pickle')
        self.data.append(pickle.load(f))
        f.close()
        
        
        ### ACA ####
        f = open(directory+'Resolution-ACA_12_0.pickle')
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
                    
                    ## We restrict only the first 4 configurations with ACA, but two first are much better
                    
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
        
        
        strOut = "### arrayResolution \n"
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
                                                

        notFound = False
        maxTry = 100
        nTry = 1
        deltaFudgeFactor = 0.05
        fudgeFactor = 1.0
        
        res , TP = self.find_array2(verbose = True)
            
        if  res[0] == '' and res[1] == '' and res[2] == '' and res[3] == '' and res[4] == '' and res[5] == '':
            notFound = True
            
        while (notFound and nTry < maxTry):
            
            nTry += 1
            notFound = False
            fudgeFactor += deltaFudgeFactor
            self.resObs *= fudgeFactor
            self.lasObs /= fudgeFactor
            res , TP = self.find_array2()
                   
            if  res[0] == '' and res[1] == '' and res[2] == '' and res[3] == '' and res[4] == '' and res[5] == '':
                notFound = True
            
            
        if nTry > 1 :
            if self.silent:
                print("# No array configuration found, fudge factor applied (Tol = %3.0f %%)"%((fudgeFactor-1.)*100.))
        
        if notFound and nTry > 1:
            if self.silent:
                print("# No array found, even with fudge factor, problem ...")
            

        strOut = ""
        pcom = ""

        if self.silent:        
            print ""
            print("### Results")

        if len(res) == 0 :
            if self.silent:
                print ",,,,"
        else:
            for s in res:
                strOut += pcom+s
                pcom = ","
        strOut+=","
        strOut+=TP

        if self.silent:        
            print strOut    
        
        return strOut
        
#====== Standalone program =========================
if __name__=="__main__":
       
    
    arg = sys.argv

    if len(arg) < 6:
        print "Arguments missing \n"
        print "The correct syntax is:"
        print "python arrayResolution.py RES (arcsec) LAS (arcsec) FREQ (GHz) Declination (Degree) Y/N (ACA) \n"
        print "Example:"
        print "python arrayResolution.py 0.2 0. 640. -53 Y    ## if LAS = 0. assumes a point source"
       
        
    else :
        a = arrayRes(arg)
        a.run()
        
        
        
    

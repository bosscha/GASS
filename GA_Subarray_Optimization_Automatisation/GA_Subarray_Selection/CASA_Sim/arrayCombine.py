#!/usr/bin/python


"""
script to run simulations with an input file containing various array to be combined (12m, ACA, TP).
Return some plots and  logs.



HISTORY:

    2013.08.10:
        - class arrayCombine
        - fct parse to read the input file.
        - fct simulation to run the simulations
        - fct analysis to produce some plot
        
        
    2013.08.11:
    
        - add concatvis, clean, componentlist
        - change the refdate for each simul
        
    2013.08.12:
        - change the project
        - modify the cleaning
        
    2013.08.13:
        - fix deconvolve
        
        
    2013.08.26:
        - adding noise option to the map (tsys-atm)
        
        
    2013.08.27:
        - adding a random seed to the noise
        
        
    2013.10.21:
        - adding hourangle for the setup parameters (default transit)
        - adding refdate  for the setup parameters.
        
    2014.10.21:
        - fixing for CASA 4.2.2
        
        
    2014.11.19:
        - simanalyzing the noisy ms if they are present
    
        
RUN:

casapy -c arrayCombine.py -i run_combine.in
"""

__version__="0.3.3@2014.11.19"
__author__ ="ALMA: SL"


import os
import sys
import UVW
import numpy as np
from math import *
import re
import random

class arrayCombine:
    
    def __init__(self):
    
       self.project = "sim"
       self.image= False
       self.analysis = False
       self.isSkymodel = False
       self.isComponentList = False
       self.inbright = ""
       self.incell = ""
       self.center = "100GHz"
       self.inwidth = "8GHz"
       self.mapsize = ""
       self.indirection = "J2000 19h00m00 -40d00m00"
       self.direction = self.indirection

       self.imsize = 256
       self.cell = '0.5arcsec'
       self.niter = 500
       self.threshold = '0.1mJy'
       self.weighting = 'natural'
       self.mask = []
       
        
   
    def parseInputFile(self,inputFile):
       "Read the input file"
       
       print("# Parsing the input file %s"%(inputFile))
       
       self.simul = []
       
       f = open(inputFile)
       
       index = -1
       for s in f:
            data = s.split()
            
            if len(data) == 3:
            
                # parameter definition
                if data[1] == "=": 
                
                    if data[0] == "project":
                        self.project = data[2]
                
                    if data[0] == "array":
                        self.simul.append({'array':data[2]})
                        index += 1
                        ## set up
                        self.simul[index]['noise'] = ''
                        self.simul[index]['user_pwv'] = 0.0
                        self.simul[index]['seed'] = 11111
                        self.simul[index]['hourangle'] = 'transit'
                        
                        
                    if data[0] == "duration":
                        self.simul[index]['duration'] = data[2]
                        
                    if data[0] == "weight" :
                        self.simul[index]["weight"] = data[2]
                        
                    if data[0] == "obsmode" :
                        self.simul[index]["obsmode"] = data[2]
                                             
                    if data[0] == "noise" :
                        self.simul[index]["noise"] = 'tsys-atm'
                        self.simul[index]["user_pwv"] = float(data[2])
                        self.simul[index]["seed"] = random.randint(1000, 100000)
                        
                    if data[0] == "declination" :
                        self.indirection = self.declinationStr(float(data[2]))
                        self.direction = self.indirection
                               
                    if data[0] == "hourangle" :
                        self.simul[index]['hourangle'] = data[2]
                        
                        
                    if data[0] == "model":
                        self.modelFile = data[2]
                        self.parseModel(data[2])
                        
                        
                    if data[0] == "image":
                        if data[2] == "True": 
                            self.image = True
                            
                    if data[0] == "imsize":
                        self.imsize = int(data[2])
                        
                    if data[0] == "cell":
                        self.cell = data[2]  
                        
                    if data[0] == "niter":
                        self.niter = int(data[2])
                        
                    if data[0] == "threshold":
                        self.threshold = data[2]
                        
                    if data[0] == "weighting":
                        self.weighting = data[2] 
                        
                    if data[0] == "mask":
                        self.mask = data[2]            
                            
                    if data[0] == "analysis":
                        if data[2] == "True":
                            self.analysis = True
                                               
   
       f.close()
       
       
       
   
    def parseModel(self,modelFile):
       "Parse the model file. If component.cl do not use it"
       
       print("# Parsing the model file %s"%(modelFile))
       
       f = open(modelFile)
       
       for s in f:
            data = s.split()
            
            if len(data) == 3:
            
                # parameter definition
                if data[1] == "=": 
                
                    if data[0] == "skymodel":
                        self.skymodel = data[2]
                        self.isSkymodel = True
                        
                    if data[0] == "inbright":
                        self.inbright = data[2]
                        
                    if data[0] == "incell":
                        self.incell = data[2]
                                     
                    if data[0] == "incenter":
                        self.incenter = data[2]
                        
                    if data[0] == "inwidth":
                        self.inwidth = data[2]
                        
                    if data[0] == "mapsize":
                        self.mapsize = data[2]
                        
                        
                    if data[0] == "componentlist":
                        self.componentlist = data[2]
                        self.isComponentList = True
                        
                    if data[0] == "compwidth":
                        self.compwidth = data[2]
                        
                                    
       
       f.close()
   
   
    def declinationStr(self,declination):
        "return the string of the coordinates"
        
        if declination <= -10:
            decString="J2000 0h00m00s %3dd00m00s"%(declination)
        elif declination <0 and declination>-10:
            decString="J2000 0h00m00s -0%1dd00m00s"%(declination)
        elif declination>=0 and declination < 10 :
            decString="J2000 0h00m00s +0%1dd00m00s"%(declination)
        elif declination >=10:
            decString="J2000 0h00m00s +%2dd00m00s"%(declination)
            
        return(decString)
        
        
        
    def summarySimulation(self):
        "Write the summary of the simulations to be performed"
        
        
        strSummary = "# Summary of the simulations \n"
        strSummary += "# Project:  %s \n"%(self.project)
        strSummary += "# Model file: %s \n"%(self.modelFile)
        strSummary += "# \n"
        strSummary += "# Model parameters: \n"
        if self.isSkymodel:
            strSummary += "# skymodel: %s \n"%(self.skymodel)
            strSummary += "# inbright: %s \n"%(self.inbright)
            strSummary += "# incell: %s \n"%(self.incell)
            strSummary += "# incenter: %s \n"%(self.incenter)
            strSummary += "# inwidth: %s \n"%(self.inwidth)
        if self.isComponentList:
            strSummary += "# componentlist: %s \n"%(self.componentlist)
            strSummary += "# compwidth: %s \n"%(self.compwidth)
            
        strSummary += "# indirection: %s \n"%(self.indirection)
        strSummary += "# mapsize: %s \n"%(self.mapsize)
        strSummary += "# \n"
        
        for arr in self.simul:
            strSummary += "# Array: %s \n"%(arr['array'])
            strSummary += "# Duration: %s \n"%(arr['duration'])
            strSummary += "# Weight: %s \n"%(arr['weight'])
            strSummary += "# Obsmode: %s \n"%(arr['obsmode'])
            strSummary += "# Reference Date: %s \n"%(arr['refdate'])
            strSummary += "# Hour angle: %s \n"%(arr['hourangle'])
            if arr['noise'] != '':
                 strSummary += "# PWV: %3.2f mm \n "%(arr['user_pwv'])
            strSummary += "# \n"
            
        print strSummary
   
    def simulation(self,param):
        "run a simulation"
        
        
        projectName = 'sim'
      

        map = 'ALMA'
        
        default(simobserve)
        
        if param['obsmode'] == 'sd':
            print("# Add 50% for the SD mapsize and a square shape")
            
            mapsi = float(re.findall("\d+.\d+|\d",self.mapsize)[0])
            mapsizeStr = re.findall('arcmin|arcsec',self.mapsize)
            mapsi *= 1.5
            mapsize = "%3.1f%s"%(mapsi,mapsizeStr[0])
            
            map = 'square'
        else:
            mapsize = self.mapsize
            
        print mapsize
        if self.isSkymodel:
            simobserve( 
                       project = projectName,
                       skymodel = "../"+self.skymodel,
                       inbright = self.inbright,
                       incell = self.incell,
                       indirection = self.indirection,
                       incenter = self.incenter,
                       inwidth  = self.inwidth,
                       direction = self.direction, 
                       antennalist = "../"+param['array'], 
                       sdantlist = "../"+param['array'],           
                       hourangle = param['hourangle'] ,
                       totaltime = param['duration'],
                       refdate = param['refdate'],
                       mapsize = mapsize,
                       maptype = map,
                       obsmode = param['obsmode'],
                       thermalnoise = param['noise'],
                       user_pwv = param['user_pwv'],
                       seed =  param['seed'])
            
        if self.isComponentList:
            simobserve(
                       project = projectName,
                       complist = "../"+self.componentlist,
                       compwidth = self.compwidth,
                       antennalist = "../"+param['array'], 
                                           
                       hourangle = 'transit',
                       totaltime = param['duration'],
                       refdate = param['refdate'],
                       mapsize = mapsize,
                       maptype = map,
                       obsmode = param['obsmode'], 
                       thermalnoise = param['noise'],
                       user_pwv = param['user_pwv'],
                       seed =  param['seed'])
    
    
    def concatVis(self):
        "Concat the interferometric MS"
        
        print("# Concat the interferometric MS")
        vis = []
        w = []
        for arr in self.simul:
            if arr['noise'] != '':
                vis1 = 'sim/sim.%s.noisy.ms'%(arr['array'].split('.')[0])
            else:
                vis1 = 'sim/sim.%s.ms'%(arr['array'].split('.')[0])
                
            if arr['obsmode'] == 'int':
                vis.append(vis1)
                w.append(float(arr['weight']))
               
        print vis 
        concat(
               vis,
               concatvis = self.project+'-full.int.ms',
               visweightscale = w)
       
    
    
    def setRefdate(self):
        "Set a different refdate"
        
        
        ref = "2014/01/"
        
        index = 0
        date = 10
        for i in self.simul:
            self.simul[index]['refdate'] = "%s%2d"%(ref,date+index)
            print self.simul[index]['refdate']
            index += 1

        
        
    def deconvolve(self):
        "run a first clean on simulations, including the full MS"
        
        
        for arr in self.simul:
            print("# Cleaning array %s"%(arr['array']))
            
            projectName = 'sim'
            if arr['noise'] != '':
                visMs = '%s/sim.%s.noisy.ms'%(projectName,arr['array'].split('.')[0])
            else:
                visMs =  '%s/sim.%s.ms'%(projectName,arr['array'].split('.')[0])
            

            print visMs
            
            simanalyze(
                        project = projectName,
                        vis = visMs ,
                        imsize = self.imsize,
                        cell   =  self.cell,
                        niter  = self.niter,
                        threshold = self.threshold,
                        weighting = self.weighting,
                        mask = self.mask
                        )
    

            

            
        
    def analysis(self,param):
        "run a first analysis"
           
        
        pass
    
    
    
    def run(self):
        "run the simulations"
        
        index = 0
        for arg in sys.argv:
            index += 1
            if arg == "-i":            
                inputFile = sys.argv[index]
    
        
        
        self.parseInputFile(inputFile)

        # Set different refdate for each simulation
        self.setRefdate()
    
        self.summarySimulation()
    
    
        ## create subdirectory project and move to it
        
        if not os.path.isdir(self.project):  
            os.mkdir(self.project)
       
        print("# Changing directory to %s"%(self.project))
        os.chdir(self.project)
        
        print("# Performing the simulations...")
        for arr in self.simul:
            print("# Array: %s \n"%(arr['array']))
            self.simulation(arr)
        
         
        ## Concat the MS
        print("# Concat the MS.")
        self.concatVis()
        
        ## Cleaning the MS
        if self.image:
            self.deconvolve()
        

#=================== Main Program =================================
#
if __name__ == "__main__":
    
    a= arrayCombine()
    a.run()





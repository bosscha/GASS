#!/usr/bin/python
# -*- coding: latin-1 -*-

"""
Classes / script to assess a given configuration following the requirements of PMG
Command is like
Array.py -f"freq (GHz)", -i "inputfile", -o "outputfile" -h (hour angle: default is -0.5 hr), -l (duration: default 1hr), -d (declination)
-h can be specified one number, but assume 1hr obs. For example, specify -2 means -2 -> -1 hr. If -l specifies, observations last that time.

HISTORY:

    2013.04.12
        - first shot
        - parsing arguments 
        
    2013.04.16:
        - The CASA part needs to be decoupled to be able to use options.
        - The simulations are tested in CASA 3.4.0
        
    2013.04.19:
        - update the casapy call
    
    2013.04.20:
        - adding parameters in the param.sim file
        
    2013.05.07:
        - add an option to allow a casa format configuration file
        - print the version
        
    2013.05.09:
        - add an option to compare input configuration with standard ones (--standard) of cycle1
        
        
    2013.07.03:
        - add an option to specify the size (arcsec) of the disk component.
        
        
    2013.09.12:
        - add the weighting for the imaging
        
        
    2013.09.24:    
        - can be used with casapy as follows : casapy -c arrayConfiguration -i input.cfg ...
        
    2014.01.02:
        - add the "-p --online" option to be used by the array configurator.
        
    
    2014.09.30:
        - Change option -c to the cycle (Default now 3)
        
    
    2014.10.08:
        - change option -c (cycle) to -x to avoid confusion with casapy -c script
        
        
    2014.11.07:
        - fix the -x option display
         
    2014.11.24:
        - add the option -t --timeRation for the computation of the time ratio if combined with standard array

    2014.12.18:
        - minor fix
        
    2014.12.23 :
        - update
    
        
    2015.03.22:
        - add the option -c N to list the N most critical pads
        
    2015.07.09:
        - change the default number of critical pads to 10.
        
    2015.07.22:
        - stable version : sync with arrayConfigurationSim
        
    2016.07.08:
        - add the option -j file.json to output a json file.
        
    2016.07.12:
        - change casapy to casa
        
    2016.08.08:
        - add the option -a (--aos) to request the current configuration at the AOS and to compute its properties.
        - change option '-s yes' to '-s'
        - add the option --rA 'List of antennas' to remove some antennas from the configuration file.
        
    2016.08.10
        - add option --twelve to select all 12-m antennas from the AOS 
        - add option --seven to select all 7-m antennas from the AOS
        
        
    2017.04.25:
        - add the option --casa to set the complete path of the casa executable
        
        
    2017.05.11:
        -updating the AOS url to https://asa.alma.cl/dashboard2-backend/service/api/antennas
        - add an option for the casa project name (--project)
        - changing the interaction with casa  to customized  param file.
        
    2017.05.12:
        - add option --cleandisk to remove some directories
        
    2017.05.14:
        - still minor bug: configuration file and project name  cannot  be  the same
        
        
    2019.01.17:
        -PRIVATE fork
        - comments go to git only
        
RUN:

python arrayConfiguration.py -h  (for help)
casapy -c arrayConfiguration.py  [option list]

For the AOS configuration and  removing the antennas  DA51, DV2 and DV3
python arrayConfiguration.py -y --rA 'DA51 DV2 DV3'

For the AOS configuration and selecting CM antennas except CM03 :
python arrayConfiguration.py -y --seven --rA 'CM03'


NOTE that for the aos only the DV and DA antennas are included (NOT PM neither CM)
"""

__version__="2.0.0@2019.01.17"
__author__ ="ALMA: SL"

import sys
sys.path.append('./')

from optparse import OptionParser
import arrayConfigurationTools as aC
import os
from time import gmtime, strftime
import json, requests
        
PAD_CFG_FILE = 'Pads.cfg'

class report:
    "Main class to produce the report"
    
    
    def __init__(self):
        
        pass
    
    
    def run(self):
        "Produce the main report"
        
        
        ## Reading the options
        
        parser = OptionParser()

        parser.add_option("-f", "--frequency", dest = "frequency", default = '100.',
                  help="Observation frequency (GHz) ")
                  
        parser.add_option("-i", "--inputfile", dest = "inputfile", 
                  help="Input file with the pad list")
        
        parser.add_option("-t", "--type", dest = "type", default = 'casa' ,
                  help="Input file type: ascii, casa")
        
                
        parser.add_option("-o","--outputfile", dest = "outputfile", default = 'arrayConfigurationResults.txt',
                          help = "Outputfile with the results")
        
        parser.add_option("-a","--hourAngle", dest = "hourangle", default = 'transit',
                          help = "Hour angle of the observation  (transit is the default)")
        
        parser.add_option("-l","--length", dest = "length", default = "1h",
                          help = "Track duration (1h is the default)")
        
        
        parser.add_option("-d","--declination", dest = "declination", default = "-23",
                          help = "Declination (-23 is the default)")
        
        parser.add_option("-k","--disk", dest = "disk", default = "1",
                          help = "Size in arcsec of the disk component")
        
        parser.add_option("-w","--weighting", dest = "weight", default = "briggs",
                          help = "Weighting (uniform, briggs, natural) to compute the array properties")
        
        parser.add_option("-p","--plot", dest = "online", default = "no",
                          help = "to produce automatic/standard plots for the online version (yes/no)")
        
        
        parser.add_option("-x","--cycle", dest = "cycle", default = "3",
                          help = "Configuration standard for Cycle nn")
        
        parser.add_option("-r","--timeRatio", dest = "timeRatio", default = "no",
                          help = "Compute the time ration if combined with standard configurations (yes|no)")
        
        
        parser.add_option("-c","--critical", dest = "critical", default = "10",
                          help = "Number of critical pads to be computed")       
        
        parser.add_option("-j","--jason", dest = "jasonFile", default = "report.json",
                          help = "Add a json file for the output.")   
        
        parser.add_option("-y","--aos", dest = "aos", default = False, action = "store_true", help = "Use the current aos configuration.")   

        
        parser.add_option("-s", "--standard", dest = "standard", default = False, action = "store_true" , help="Compare with a standard configuration")

        
        parser.add_option("--rA", dest = "removeAntenna", default = "", help = "Only for AOS ! Remove antenna from the configuration file. Enter  the list as a string") 

        
        parser.add_option("--twelve", dest = "all12m", default = False, action = "store_true" , help = "Only for AOS ! Select all 12-m antennas (DV, DA, PM)") 

        parser.add_option("--seven", dest = "all7m", default = False, action = "store_true" , help = "Only for AOS ! Select all 7-m antennas (CM)") 

        parser.add_option("--casa", dest = "casaExec", default = "casa", help = "Executable for casa (/path/casa) ") 
        
        parser.add_option("--project", dest = "project", default = "sim", help = "project name for the CASA simulation") 


        parser.add_option("--cleandisk", dest = "cleanDir", default = False, action = "store_true" , help = "To remove some directories and reduce  disk space")

        (optarr, args) = parser.parse_args()
        
        
        if  optarr.standard :
            optarr.standard = 'yes'
        else :
            optarr.standard = 'no'
                   
        if   optarr.aos :
            print("## Obtaining the AOS configuration ...")
            aosCurrent = aosConfiguration()
            
            
            if optarr.removeAntenna != "":
                antennaToRemove = optarr.removeAntenna.split()
                print("#### Do not consider : %s"%antennaToRemove)
                optarr.inputfile = aosCurrent.getAOS(antennaToRemove, all12 = optarr.all12m, all7 = optarr.all7m)               
            else :
                optarr.inputfile = aosCurrent.getAOS(all12 = optarr.all12m, all7 = optarr.all7m)

            optarr.type = 'casa'
                      
        ## cleaning disk space option   
        if optarr.cleanDir :
            cleanDir = 'True'
        else:
            cleanDir = 'False'
            
            
        print optarr.standard
        
        print("## Simulation of an array configuration")
        print("## version: %s"%(__version__))
        print("## Input file: %s"%(optarr.inputfile))
        print("## Output file: %s"%(optarr.outputfile))
        print("## Frequency (GHz): %s"%(optarr.frequency))
        print("## Declination (Degree): %s"%(optarr.declination))
        print("## Observation Start (Hour Angle): %s"%(optarr.hourangle))
        print("## Duration: %s"%(optarr.length))
        print("## Standard: %s"%(optarr.standard))
        print("## Time Ratio: %s"%(optarr.timeRatio))
        print("## Disk component: %s (arcsec)"%(optarr.disk))
        print("## Weighting: %s "%(optarr.weight))
        print("## Online: %s "%(optarr.online))
        print("## Critical Pads: %s "%(optarr.critical))
        print("## Cycle: %s"%(optarr.cycle))
        print("## Project: %s"%(optarr.project))
        print("## Jason File:  %s"%(optarr.jasonFile))
        print("## Clean Directories:  %s"%(cleanDir))       
        print("##")
    
        ## Creating the CASA configuration file
        a = aC.ArrayConfigurationCasaFile(padPositionFile = PAD_CFG_FILE)
        
        
        listpad = []
        listant = []
        
        
        f = open(optarr.inputfile)
        
        for pad in f:
            if pad[0] != "#":
                data = pad.split()
                if len(data) == 2:
                    listant.append(data[0])
                    listpad.append(data[1])
            
            
        f.close()
        
        Nant = len(listpad)
        
        # print listpad
        if optarr.type == 'ascii':
            a.createCasaConfig(optarr.inputfile,listPads = listpad)
            casaFile = optarr.inputfile+'.cfg'
            
        elif optarr.type == 'casa':
            casaFile = optarr.inputfile
                  
        
        # Info 
        
        ai = aC.ArrayInfo(casaFile)
        ai.stats()
         
        
        # produce input for CASA scripts
        simParamFile = "%s.par"%(optarr.inputfile)
        f = open(simParamFile,"w")
        
        f.write("file         = %s \n"%(casaFile))
        f.write("inputfile    = %s \n"%(optarr.inputfile))
        f.write("reportFile   = %s \n"%(optarr.outputfile))
        f.write("Nant         = %s \n"%(Nant))
        f.write("frequency    = %s \n"%(optarr.frequency))
        f.write("declination  = %s \n"%(optarr.declination))
        f.write("hourAngle    = %s \n"%(optarr.hourangle))
        f.write("duration     = %s \n"%(optarr.length))
        f.write("minBl        = %f \n"%(ai.minBl))
        f.write("maxBl        = %f \n"%(ai.maxBl))
        f.write("rms          = %f \n"%(ai.rms))
        f.write("standard     = %s \n"%(optarr.standard))
        f.write("timeRatio    = %s \n"%(optarr.timeRatio))
        f.write("disk         = %s \n"%(optarr.disk))
        f.write("weighting    = %s \n"%(optarr.weight))
        f.write("cycle        = %s \n"%(optarr.cycle))
        f.write("online       = %s \n"%(optarr.online))
        f.write("criticalPads = %s \n"%(optarr.critical))
        f.write("jasonFile    = %s \n"%(optarr.jasonFile))
        f.write("projectName  = %s \n"%(optarr.project))
        f.write("cleanDir     = %s"%(cleanDir)) 
        
        f.close()
        
        print("## Simulation parameters written")

        ## run the simulations
        
        cmd = "casa  --nogui -c arrayConfigurationSim.py -i %s"%(simParamFile)
        cmdalt = "casapy  --nogui -c arrayConfigurationSim.py -i %s"%(simParamFile)
        
        ## change the casa executable if needed
        if optarr.casaExec != "casa":
            cmd = "%s  --nogui -c arrayConfigurationSim.py -i %s"%(optarr.casaExec, simParamFile)
        
        try:
            os.system(cmd)
        except:
            os.system(cmdalt)


class aosConfiguration:
    "class to obtain and to treat the AOS configuration"
    
    def __init__(self):
        
        self.aosFile = "aos-%s"%(strftime("%Y%m%d-%Hh%Mm%Ss", gmtime()))

        
    def getAOS(self,  antennaToRemove = [] , all12 = False, all7 = False):
        "Get the current AOS configuration (only 12-m yet)"
        
        try:
             resp = requests.get('https://asa.alma.cl/dashboard2-backend/service/api/antennas')

        except:
            print("### Error on getting the current AOS configuration ....")
            
        
        dataAOS = json.loads(resp.content)   
        
        
        faos = open(self.aosFile,'w')
        
        listPossibleAnt = ['DV', 'DA']
        
        if all12 :
            listPossibleAnt.append('PM')
            
        elif all7 :
            listPossibleAnt = ['CM']
        
        
        for key in dataAOS:
            if key['antenna-functional-status'] == 'AVAILABLE' and (key['name'][0:2] in listPossibleAnt)  \
              and key['name'] not in antennaToRemove:

                faos.write("%s  %s \n"%(key['pad'],key['name']))
        
        
        faos.close()
        
        a = aC.ArrayConfigurationCasaFile(padPositionFile = PAD_CFG_FILE)
        a.createCasaConfig(self.aosFile)
        
        return(self.aosFile+'.cfg')
        
#=================== Main Program =================================
#
if __name__ == "__main__":
    
    a = report()
    a.run()
    print "### Done"
    
    
    

#!/usr/bin/python


"""
Set of methods to create an JSON input file for ACDC.

HISTORY:

    2014.10.01:
        - first shot

    2014.10.21:
        - add rms keyword (-b)
        
    2014.10.30:
        - add result file (-f)
    

RUN

"""

 
__authors__="ALMA: SL"
__version__="0.2.0@2014.10.21"

from optparse import OptionParser
import json

ALMA_LATITUDE = -23.0

class acdc2json:
    "Class to create the input json file for ACDC"
    
    def __init__(self, file = "input.json"):

        self.fileJson = file
        self.data = {}
        
        
    def set_version(self,version):
        "Set the version"
        
        self.data['version'] = {'value': version,'doc': 'Version of the input file'}
        
        
    def set_nadd(self,nadd):
        "Number of antennas to add to the initial configuration"
        
        self.data['nadd'] = {'value': nadd,'doc': 'number of antennas to put on the destination pads'}
        
        print self.data
        
       
    def set_nremove(self,nremove):
        "Number of antennas to remove  from the initial configuration"
        
        self.data['nremove'] = {'value': nremove,'doc': 'number of antennas to remove from initial config'}     


    def set_sdec(self,sdec):
        "Source declination in degrees"
        
        self.data['sdec'] = {'value': sdec,'doc': 'Source declination in degree'}     


    def set_resol(self,resol):
        "resolution in arcsec at 345GHz, , range=[0.01, 2.0]"
        
        self.data['resol'] = {'value': resol,'doc': 'resolution in arcsec at 345GHz, , range=[0.01, 2.0]'}    
        
    
    def set_rms(self, rms):
        "RMS 1D to filter the best config"
        
        self.data['baseline_RMS'] = {'value': rms ,'doc': 'RMS of the baseline lengths wanted for the configuration in meters, range=[50.0, 10000.0]'}   
        

    def set_resol_tol(self,resol_tol):
        "resolution tolerance around the specification in percent"
        
        self.data['resol_tol'] = {'value': resol_tol,'doc': 'resolution tolerance around the specification in percent'}     


    def set_elev_tol(self,elev_tol):
        "elevation tolerance allowed around the specification in degrees "
        
        self.data['elev_tol'] = {'value': elev_tol,'doc': 'elevation tolerance allowed around the specification in degrees '}     


    def set_sitelat(self,sitelat):
        "site latitude in degrees, range=[-90, 90]"
        
        self.data['sitelat'] = {'value': sitelat,'doc': 'site latitude in degrees, range=[-90, 90]'}     


    def set_ncells(self,ncells):
        " number of cells along one axis in the grid (the total is ncells^2"
        
        self.data['ncells'] = {'value': ncells,'doc': 'number of cells along one axis in the grid (the total is ncells^2 '}     


    def set_nwrite(self,nwrite):
        "number of configs to write in the best list file"
        
        self.data['nwrite'] = {'value': nwrite,'doc': 'number of configs to write in the best list file'}     

    def set_resultfile(self,resultfile):
        "name of the file in which to write the best configs"
        
        self.data['resultfile'] = {'value': resultfile,'doc': 'name of the file in which to write the best configs'}     
   
    
    def statusPads(self, padname):
        "find the status of the antenna on the padname"
        
        fants = open(self.ants)
        
        status = 'free'
        
        for line in fants:
            dat = line.split()
            
            if dat[0] == padname :
                status = dat[1]
                
        fants.close
        
        return(status)
            
        
    def set_pads(self):
       "Set the configuration of Pads for ACDC"
       
       self.data['pads'] = {"doc" : "Here one may put some informations relative the specific set of pads which is presented to ACDC",
                            "value" : {"doc" : {"name" : "name of the pad",
                                                "position" : "coordinate of the pads (X, Y) with respect to any reference center. Must be numeric real values",
                                                "status" : "qualifies the status this pad and consequently the action which can be done on it. Must be one string in {'free', 'occupied', 'keep_occupied'}"
          }, "list" : [] }}
       
 
       
       fpads = open(self.pads)

 
       for line in fpads:
           if line[0] != '#':
               
               dat = line.split()
               x = float(dat[0])
               y = float(dat[1])
               z = float(dat[2])
               padName = dat[4]
               status = self.statusPads(padName)
               
               padDict = {"name" : padName, "position" : {"x": x, "y": y}, "status" : status }
               
               self.data['pads']['value']['list'].append(padDict)
       
       fpads.close()
       
       
       
 
 
 
        
    def run(self):
        "run the script with the options"
        
        
        parser = OptionParser()

        parser.add_option("-p", "--pads", dest = "pads", default = 'pads.cfg',
                  help="Pads.cfg file to target (CFG format) ")
        
        parser.add_option("-a", "--ant", dest = "ants", default = 'none',
                  help="pads with antennas (ascii file format PAdS OPTION")
        
        parser.add_option("-o", "--output", dest = "outputfile", default = 'input.json',
                  help="Output Json file")
        
        parser.add_option("-n", "--nadd", dest = "nadd", default = "0",
                  help="Number of antenna to add")
                
        parser.add_option("-r", "--nremove", dest = "nremove", default = "0",
                  help="Number of antenna to remove")
        
        parser.add_option("-x", "--resolution", dest = "resol", default = "-1",
                  help="Spatial resolution at 345 GHz in arcsec") 
        
        parser.add_option("-b", "--blrms", dest = "rms", default = "-1",
                  help="RMS  of the baseline length") 
        
        
        parser.add_option("-f", "--file", dest = "resultfile", default = "bestconfigs.txt",
                  help="result file name (default bestconfigs.txt)") 
        
        
        
                
        (optarr, args) = parser.parse_args()

        
        print("## ACDC input file (Json) of an array configuration")
        print("## version: %s"%(__version__))
        print("##")        
        print("## Allowed Pad file: %s"%(optarr.pads))
        print("## Antenna file: %s"%(optarr.ants))
        print("## Number of antenas to add: %s"%(optarr.nadd))
        print("## Number of antenas to remove: %s"%(optarr.nremove))
        print("## RES at 345 GHz : %s arcsec"%(optarr.resol))
        print("## RMS of the BL length : %s m"%(optarr.rms))
        print("## Output file: %s"%(optarr.outputfile))
        print("## Result file: %s"%(optarr.resultfile))
        print("##")
        
        self.pads = optarr.pads
        self.ants = optarr.ants
        self.fileJson = optarr.outputfile
        
        self.set_version(1)
        self.set_nadd(int(optarr.nadd))
        self.set_nremove(int(optarr.nremove))     
        self.set_sdec(-53.)
        if optarr.resol != "-1":
            self.set_resol(float(optarr.resol))
        if optarr.rms != "-1":
            self.set_rms(float(optarr.rms))
        self.set_resol_tol(15.0)
        self.set_elev_tol(10.0)
        self.set_sitelat(ALMA_LATITUDE)
        self.set_ncells(5)
        self.set_nwrite(100)
        self.set_resultfile(optarr.resultfile)
        self.set_pads()
        
        
          
        
        
        # print json.dumps(self.data, skipkeys=True, indent = 2, separators=(',', ': '))
        
        fout = open(optarr.outputfile, "w")
           
        fout.write(json.dumps(self.data, skipkeys=True, indent = 2, separators=(',', ': ')))
        
        fout.close()
        
        
        
        
#=================== Main Program =================================
#
if __name__ == "__main__":
    
    a = acdc2json()
    a.run()
    print "### Done"
        
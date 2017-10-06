"""
Script to run on the virtual machine in order to get the array information for a given uid (EB)


HISTORY:
    2016.10.24:
        - first shot from the Alvaro's example
        
    
    2016.10.25:
        - adding the parsing of xml files
        - creating the configuration file for the array
        
    2017.03.30:
        - add asdmAnalysis class to provide some method for the analysis
        
    2017.03.31:
        - updating the findHAminmax
        - copy and adapt computeLST from analysisUtils
        
        
    2017.04.03:
        - bug fixing...
        
        
    2017.04.04:
        - bug fixing
        - add a run method.
        - add the frequency dependence.
    
        
    2017.04.05:
        - add the case of ACA antennas...
        
    
    2017.05.11:
        - merging version with service version (AA)
        
    2017.05.12:
        - clean some directories.
        - add a specific project for each run
        
    2017.05.14:
        - debugging cfg file creating

    2017.06.06:
        - add a method to discard flagged data/antennas.
        - start coding to use direct position of the antennas (w/o Pads.cfg file)
        
        
    2017.07.07:
        - fixing bugs
        - if directionSource = [0.0, 0.0] it is probably a solar system object. As a temporary fix we set HA to 0.
        
"""

__version__="0.4.3@2017.07.07"
__author__ ="ALMA: SL, AA"


from lxml import etree
from lxml.etree import XMLParser
from asdf.archive import asdm_export

import numpy as  np
import arrayConfigurationTools as act
import os

from optparse import OptionParser

#aaguirre
from wsgiref.simple_server import make_server
from pyramid.config import Configurator
from pyramid.response import Response
from pyramid.response import FileResponse

from asdf.utils import is_valid_uid
from pyramid.httpexceptions import HTTPBadRequest

p = XMLParser(huge_tree=True)

ALMA_LONGITUDE=-67.754748 # ICT-4143,  -67.754694=JPL-Horizons,  -67.754929=CASA observatories
ALMA_LATITUDE=-23.029211  # ICT-4143,  -23.029167=JPL-Horizons,  -23.022886=CASA observatories
BB_NAME = ['BB_1', 'BB_2','BB_3','BB_4']


def readOptions():
    "Read the CL options ..."

    parser = OptionParser()

    parser.add_option("-u", "--uid", dest = "uid", default = 'none',
                  help="uid of the EB. Ex. : uid://A002/Xbd3836/X5416")
    
    
    
    (optarr, args) = parser.parse_args()            
    
    options = {}    
    options['uid']  = optarr.uid

    return(options)


class asdmAnalysis:
    "class to provide some tools for the analysis of the metadat in the asdm"
    
    def __init__(self, uid):
        
        data = asdm_export(uid,['Antenna', 'Station', 'SpectralWindow', 'Scan', 'Source','Flag'])
    
        self.antenna_xml        = data.get('Antenna').xml
        self.station_xml        = data.get('Station').xml
        self.spw_xml            = data.get('SpectralWindow').xml
        self.scan_xml           = data.get('Scan').xml
        self.source_xml         = data.get('Source').xml
        self.spectralWindow_xml = data.get('SpectralWindow').xml
        self.flag_xml           = data.get('Flag').xml

        # aaguirre
        self.conf_filename = uid.replace(":","_").replace("/","_")
        
    
    def createConfigurationFile(self):
        "create the uid.cfg configuration file"
    
        dish = []
        antName = []
        stationId = []
        padName = []
    
        xml_object = etree.fromstring(self.antenna_xml)
        rows_ant = xml_object.findall('row')
     
    ## testing ...
        ftest = open("test.cin","w")
    
        for row in rows_ant :
            res = row.find('dishDiameter').text.strip()
            dish.append(res)
        
            res = row.find('name').text.strip()
            antName.append(res)
        
            res = row.find('stationId').text.strip()
            
            position = row.find('position').text.strip()
            ftest.write(position)
            
            stationId.append(res)
        
            padName.append('none')
        
        ##testing
        ftest.close()
        
        ## find the match of stationId with Padname
        xml_object = etree.fromstring(self.station_xml)
        rows_station = xml_object.findall('row')   
    
        for row in rows_station:
            pad  = row.find('name').text.strip()
            id = row.find('stationId').text.strip()
        
            for i in range(len(dish)):
                if stationId[i] == id :
                    padName[i] = pad
                            
    
        ## creation of the array configuration file
        aCasa = act.ArrayConfigurationCasaFile()
        listpad = []
    
    
        ## for 12-m array
        for i in range(len(dish)):
            if int(float(dish[i])) == 12 :
                listpad.append(padName[i])
                
            if int(float(dish[i])) == 7 :
                listpad.append(padName[i])
            
        # rm any possible uid.cfg

        # aaguirre
        cmd = "rm %s.cfg" % self.conf_filename
        os.system(cmd)
        
        aCasa.createCasaConfig(self.conf_filename, listPads = listpad)
        

        
        return(0)
    

    def findFlaggedAntennas(self, startTime, endTime, fractionTime = 0.5):
        """
        Find the antennas to be flagged. 
            - startTime, endTime: time window to analyze (ALMA time)
            - fractionTime: minimum fraction time to consider an antenna unflagged
        """
        
        ant_flag = {}
        
        xml_object = etree.fromstring(self.flag_xml)
        rows_flag = xml_object.findall('row')
        
        
        for row in row_flag:
            start = row.find('startTime').text.strip()
            end   = row.find('endTime').text.strip()
            antId = row.find('antennaId').text.strip()
    
    
    def ComputeLST(self, mjdsec, longitude = ALMA_LONGITUDE):
        """
        Computes the LST (in hours) for a specified time and longitude. 
        The input longitude is in degrees, where east of Greenwich is > 0.

        mjdsec: MJD seconds

        """
    
    
        JD =  mjdsec / 86400. + 2400000.5       
        T = (JD - 2451545.0) / 36525.0
        sidereal = 280.46061837 + 360.98564736629*(JD - 2451545.0) + 0.000387933*T*T - T*T*T/38710000.

        # now we have LST in Greenwich, need to scale back to site
            
        sidereal += longitude
        sidereal /= 360.
        sidereal -= np.floor(sidereal)
        sidereal *= 24.0
    
        if (sidereal < 0):
            sidereal += 24
        if (sidereal >= 24):
            sidereal -= 24

        return(sidereal)

   
    def findSourceDirectionSpw(self,source):
        "find the direction of a source"
        
        xml_object = etree.fromstring(self.source_xml)
        rows_source = xml_object.findall('row')
        
        found = False
        directionSource = [0.,0.]
        spwIdList = []
        
        for row in rows_source:
            sourceName = row.find('sourceName').text.strip()           
            direction  = row.find('direction').text.strip()
            spwId      = row.find('spectralWindowId').text.strip()
            
            if sourceName == source:
                found = True
                directionSource = direction.split(" ")[2:4]
                spwIdList.append(spwId)
                
            
        print spwIdList
                
        
        return(found, directionSource, spwIdList)
    
    
    def findFrequency(self, spwId):
        "Find the frequency of the source using the spwId"
        
        refFreqMean = 0.
        nFreq       = 0
        
        xml_object = etree.fromstring(self.spectralWindow_xml)
        rows_spw   = xml_object.findall('row')
        
        for row in rows_spw:
            spwIdRow   = row.find('spectralWindowId').text.strip()
            refFreq    = float(row.find('refFreq').text.strip()) / 1e9
            bbName     = row.find('basebandName').text.strip()
            
            if (spwIdRow in spwId) and (bbName in BB_NAME):
                refFreqMean += refFreq
                nFreq += 1.
                
        refFreqMean /= nFreq
        
        print("frequency")
        print refFreqMean
        
        return(refFreqMean)
    
    
   
    def findHAminmax(self):
        "Find the HA min and max on the first science target ..."
        
        
        xml_object = etree.fromstring(self.scan_xml)
        rows_scan = xml_object.findall('row')
        
        
        targetSource = ''
        foundSource = False
        
        for row in rows_scan:
            intent    = row.find('scanIntent').text.strip().split(' ')           
            startTime = float(row.find('startTime').text.strip()) / 1e9
            endTime   = float(row.find('endTime').text.strip())   / 1e9
            source    = row.find('sourceName').text.strip()

            print "Intent"
            print intent
            
            print("Startime")
            print startTime
            
            
            print("Endtime")
            print endTime
            
            print("source")
            print source
            
            
            if 'OBSERVE_TARGET' in intent:
                if not foundSource:
                    targetSource = source
                    foundSource = True
                    
                    time0 = startTime
                    time1 = endTime
                    
                if foundSource and source == targetSource:
                    if startTime < time0:
                        time0 = startTime
                        
                    
                    if endTime > time1:
                        time1 = endTime
        
        print("Target source")                
        print targetSource
        print time0
        print time1
                        
        foundDirection , directionSource , spwIdList = self.findSourceDirectionSpw(targetSource)
        
        if not foundDirection :
            return(-99,-99, -99, -99)
        
        ## find the frequency of the source
        frequency = self.findFrequency(spwIdList)
        
    
        ## from analysisUtils
        minmjd = time0 / 86400.
        maxmjd = time1 / 86400.
        
        
        minLST = self.ComputeLST(minmjd * 86400.)
        maxLST = self.ComputeLST(maxmjd * 86400.)
        
        print("minlst")
        print minLST
        print directionSource 
        
        minHourAngle = minLST - float(directionSource[0]) * 12 / np.pi
        maxHourAngle = maxLST - float(directionSource[0]) * 12 / np.pi
        
        if(minHourAngle > 12.):
            minHourAngle -= 24.
            maxHourAngle -= 24.
        
        ## fix for solar system object but with wrong HA
        if  float(directionSource[0]) == 0. and  float(directionSource[1]) == 0.:
            a0 = (maxHourAngle - minHourAngle) / 2.
            minHourAngle = -a0
            maxHourAngle =  a0
            
        
        return(minHourAngle , maxHourAngle, float(directionSource[1])* 180. / np.pi, frequency)
    
    
    def run(self):
        "run the analysis on a uid (metadata"
        
            
        ## create the configuration file ...   
        print("## Creating array configuration file...")
        self.createConfigurationFile()
        
        ### find HAmin and HAmax for the first science target ..
        haMin , haMax , declinationSource, frequencySource = self.findHAminmax()

        print("##")
        print("# hamin, hamax, declination")
        print haMin
        print haMax
        print declinationSource
        
        
        # run the simulation and create the json file uid.json. Compute teh full list of critical antennas.
    
        ha = (haMin + haMax) / 2.
        duration = haMax - haMin
        ## aaguirre
        
        cmd = "python arrayConfiguration.py -i %s.cfg -t casa -c 100 -f %f -d %f -a %fh -l %fh --j %s.json --project %sSim  --cleandisk  "%(self.conf_filename,frequencySource, declinationSource, ha, duration, self.conf_filename, self.conf_filename)  
        print cmd
        os.system(cmd)



def array_eb(request):
    uid = request.params.get('uid',None)
    if uid is None: raise HTTPBadRequest('No valid UID provided')
    
    if not is_valid_uid(uid): raise HTTPBadRequest('No valid UID provided')
    aA = asdmAnalysis(uid)
    aA.run()
    path = '/home/sleon/testing/' +  uid.replace(":","_").replace("/","_") + '.json'
    response = FileResponse(
          path,
          request=request,
          content_type='application/json'
        )

    return response

#######################################################################################################   
#######################################################################################################    
if __name__ == '__main__':
    # uid = 'uid://A002/Xb25e1a/X166e1'
    # uid = 'uid://A002/Xbd3836/X5416'
    # uid = 'uid://A002/Xbebcb7/X1c8b'   ## ACA

    
    ## options..
    opt = readOptions()
    uid  = opt['uid']
    
    if uid == 'none':
        print("## No uid...")
        print("## Done")
    
    else:
        print uid
        aA = asdmAnalysis(uid)
        aA.run()
        print("### Done")
        
        
    """
    config = Configurator()
    config.add_route('array_eb', '/array-eb')
    config.add_view(array_eb, route_name='array_eb')
    app = config.make_wsgi_app()
    server = make_server('0.0.0.0', 8080, app)
    server.serve_forever()
    """

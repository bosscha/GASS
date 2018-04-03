#!/usr/bin/python


"""
class/method to assess the imaging quality of an array configuration. It is complementary to the metrics produced by the 
arrayConfiguration.py script.


HISTORY:
    2015.09.08:
        - first shot with only the comparing part without including the simulations.
    
    2015.09.10:
        - implement first metric [(o-i) / noise]
        - !! Still missing a proper comparison with different spatial resolution
        - regridding output image to the input image in temp.regridded.image

RUN:
> casapy -c arrayConfigurationImaging.py [OPTIONS]
> casapy -c arrayConfigurationImaging.py -i image1 -o image2 -s 3.0 -n 5e-3s


"""


__version__="0.0.2@2015.09.10"
__author__ ="ALMA: SL"

from optparse import OptionParser
import numpy as np
import math

RAD2DEG = 180./ math.pi

class fidelity:
    
    def __init__(self):
        
        pass
    
    def run(self):
        "parse the option and run the computation"
        
        
        ## Reading the options
        
        parser = OptionParser()

        parser.add_option("-i", "--imagein", dest = "imageIn", default = 'imageIn',
                  help="image IN sky) in the CASA format ")
                  
        parser.add_option("-o", "--imageout", dest = "imageOut",default='imageOut', 
                  help="Input file with the pad list")
        
        parser.add_option("-s", "--sigma", dest = "sigma",default= '5', 
                  help="minimum sigma to analyze the imageOut")
        
        parser.add_option("-n", "--noise", dest = "noise",default= '-1', 
                  help="Set the noise of the imageOut. If not set, it will estimate it")
        
        parser.add_option("-r", "--rnax", dest = "rmax",default= '10', 
                  help="Set the maximum length from the center to consider the two images")
        
        parser.add_option("-c", "--fidelity", dest = "script", default = 'arrayConfigurationImaging.py',
                  help="casa script for the fidelity == arrayConfigurationImaging.py")
        
        parser.add_option("-k", "--knorm", dest = "knorm",default= 'no' , 
                  help="Normalisation of the two images to the central part value (yes|no)")
        
        parser.add_option("-g", "--gridding", dest = "gridding",default= 'yes' , 
                  help="Regridding the output image to the coord. of the input (yes|no)")
        
        (optarr, args) = parser.parse_args()
        
        
        ## processing the options
        ###
        
        print("### In Testing !!")
        imageIn  = optarr.imageIn
        imageOut = optarr.imageOut
        
        noise = float(optarr.noise)
        sigma = float(optarr.sigma)
        rmax  = float(optarr.rmax)
        
        if optarr.knorm == 'no' :
            normalisation = False
        else:
            normalisation = True
        
        if optarr.gridding == 'yes' :
            gridding = True
            self.griddingInput(imageIn, imageOut)
        else:
            gridding = False
        
        if noise == -1.0:
            noise = self.estimateNoise(imageOut)
        
        
        
        self.imagingMetric1(imageIn, imageOut, rmax, gridding)
        
        
    def griddingInput(self,imageIn, imageOut):
        "We regrid the input on the output, the result goes in temp.regridded.image"
        
        imres = "temp.regridded.image"
        
        ia.open(imageOut)
        mycs = ia.coordsys()
        ia.close()
        
        ia.open(imageIn)  
        imrr = ia.regrid(outfile=imres, csys=mycs.torecord(),  
                  shape=ia.shape(), overwrite=true)
          
        mycs.done()  
        imrr.done()  
        ia.close()
        
        
       
    def estimateNoise(self,image):
        "Estimate the noise in an image"
        
        return(0.)
    
       
        
    def imagingMetric1(self,image1, image2, rmax, gridding):
        """
        Compare the image2 produced by a transform of image1. The sigma set the rms level to consider on the image2
        Return a distance between image2 and image1.
        Does not check if the pixel size are identical...
        """
        
        if gridding:
            image1 = 'temp.regridded.image'
            
        im1h = imhead(imagename= image1,mode="list",hdkey="",hdvalue="",verbose=False)
        im2h = imhead(imagename= image2,mode="list",hdkey="",hdvalue="",verbose=False)
        
        print("### Compute metric1")
        print("## Image in    : %s"%(image1))
        print("## Image out   : %s"%(image2))
        print("## Rmax (pixel): %d"%(rmax))
        
        
        xc2 = im2h['shape'][0] / 2
        yc2 = im2h['shape'][1] / 2

        print xc2, yc2

        ## We extract the same array on both image
        
        ia.open(image1)
        II = ia.getchunk([xc2-rmax,yc2-rmax],[xc2+rmax,yc2+rmax])
        im1 = np.zeros((len(II),len(II)))
        print len(II)
        
        for i in range(len(II)):
            for j in range(len(II)):
                im1[i,j] = II[i][j]
        ia.close()
        
        ia.open(image2)
        OO = ia.getchunk([xc2-rmax,yc2-rmax],[xc2+rmax,yc2+rmax])
        im2 = np.zeros((len(OO),len(OO)))
        
        print len(OO)
        
        for i in range(len(OO)):
            for j in range(len(OO)):
                im2[i,j] = OO[i][j]
        ia.close()
        
        
        im1_2   = np.sum(np.dot(im1,im1))
        im2_2   = np.sum(np.dot(im2,im2))
        im1_im2 = np.sum(np.dot(im1,im2))
            
        Q = im1_im2 / math.sqrt(im1_2 * im2_2)
        
        print im1_2 , im2_2 , im1_im2
        print Q
        diffIm = math.acos(Q) * RAD2DEG

        
        print("## Fidelity 1 [0--1]: %2.4f  (%2.2f degrees)"%(abs(Q), diffIm))
        
        return(Q)
        
        

              
              
              
    
        
#=================== Main Program =================================
#
if __name__ == "__main__":
    
    print("### run the fidelity metrics")
    a = fidelity()
    a.run()
    print "### Done"
    
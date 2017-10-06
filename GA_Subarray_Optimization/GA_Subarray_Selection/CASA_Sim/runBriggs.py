## plot the resolution for different Briggs factor
## C40-4.cfg // 100 GHz // 1h // -23d
## -2 Uniform
## 2 Natural


import numpy as np
import os
import math

import pickle


briggs = np.arange(0,101.)/25 -2.0
resolution = []

for R in briggs:
  os.system("rm -Rf toto*")
  print R
  
  clean(vis="./sim/sim.tempCasa.ms/",imagename="toto",outlierfile="",field="",spw="",selectdata=False,timerange="",uvrange="",antenna="",scan="",observation="",mode="mfs",gridmode="",wprojplanes=1,facets=1,cfcache="cfcache.dir",painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=500,gain=0.1,threshold="0.0mJy",psfmode="clark",imagermode="csclean",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,interactive=False,mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=128,cell="0.1arcsec",phasecenter="",restfreq="",stokes="I",weighting="briggs",robust=R,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,allowchunk=False)
  
  ia.open("toto.image")
  data = ia.restoringbeam()
  
  resolution.append(math.sqrt(data['major']['value']*data['minor']['value']))
  ia.close()
  
  


f = open("brigss.pickle","w")
pickle.dump([briggs,resolution],f)
f.close()

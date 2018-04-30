#!/usr/bin/python

"""
Compute and display the minor/major beam resolution for different array configuration


"""

import sys

# sys.path.insert(0,'/home/stephane/git/ALMA/ALMA/ArrayConfiguration/')
#sys.path.insert(0,'/home/stephane/workspace/AIV/science/qa2')  ## psf class


from arrayConfigurationReport import *
import pickle

b = beam()

dec0, minor0, major0, las0 = b.curveBeamResolutiongDeclination('C43-1.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-1.pickle', 'w')
data = [dec0,minor0,major0,las0]
pickle.dump(data,f)
f.close()

dec1, minor1, major1, las1 = b.curveBeamResolutiongDeclination('C43-2.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-2.pickle', 'w')
data = [dec1,minor1,major1,las1]
pickle.dump(data,f)
f.close()

dec2, minor2, major2, las2 = b.curveBeamResolutiongDeclination('C43-3.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-3.pickle', 'w')
data = [dec2,minor2,major2,las2]
pickle.dump(data,f)
f.close()

dec3, minor3, major3, las3 = b.curveBeamResolutiongDeclination('C43-4.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-4.pickle', 'w')
data = [dec3,minor3,major3,las3]
pickle.dump(data,f)
f.close()

dec4, minor4, major4, las4 = b.curveBeamResolutiongDeclination('C43-5.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-5.pickle', 'w')
data = [dec4,minor4,major4,las4]
pickle.dump(data,f)
f.close()

dec5, minor5, major5, las5 = b.curveBeamResolutiongDeclination('C43-6.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-6.pickle', 'w')
data = [dec5,minor5,major5,las5]
pickle.dump(data,f)
f.close()

dec6, minor6, major6, las6 = b.curveBeamResolutiongDeclination('C43-7.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-7.pickle', 'w')
data = [dec6,minor6,major6,las6]
pickle.dump(data,f)
f.close()

dec7, minor7, major7, las7 = b.curveBeamResolutiongDeclination('C43-8.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-8.pickle', 'w')
data = [dec7,minor7,major7,las7]
pickle.dump(data,f)
f.close()

dec8, minor8, major8, las8 = b.curveBeamResolutiongDeclination('C43-9.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-9.pickle', 'w')
data = [dec8,minor8,major8,las8]
pickle.dump(data,f)
f.close()

dec9, minor9, major9, las9 = b.curveBeamResolutiongDeclination('C43-10.cfg','2h',-80.,46.,3.)
f = open('Resolution-C43-10.pickle', 'w')
data = [dec9,minor9,major9,las9]
pickle.dump(data,f)
f.close()

dec10, minor10, major10, las10 = b.curveBeamResolutiongDeclination('ACA-std.cfg','2h',-80.,46.,3.)
f = open('Resolution-ACA-std.pickle', 'w')
data = [dec10,minor10,major10,las10]
pickle.dump(data,f)
f.close()



fig = pl.figure()
ax = fig.add_subplot(111)

ax.plot(dec0,major0,'k-')
ax.plot(dec0,minor0,'k-')
ax.fill_between(dec0,major0,minor0,facecolor='yellow')

ax.plot(dec1,major1,'k-')
ax.plot(dec1,minor1,'k-')
ax.fill_between(dec1,major1,minor1,facecolor='blue')

ax.plot(dec2,major2,'k-')
ax.plot(dec2,minor2,'k-')
ax.fill_between(dec2,major2,minor2,facecolor='green')

ax.plot(dec3,major3,'k-')
ax.plot(dec3,minor3,'k-')
ax.fill_between(dec3,major3,minor3,facecolor='red')

ax.plot(dec4,major4,'k-')
ax.plot(dec4,minor4,'k-')
ax.fill_between(dec4,major4,minor4,facecolor='yellow')

ax.plot(dec5,major5,'k-')
ax.plot(dec5,minor5,'k-')
ax.fill_between(dec5,major5,minor5,facecolor='cyan')

ax.plot(dec6,major6,'k-')
ax.plot(dec6,minor6,'k-')
ax.fill_between(dec6,major6,minor6,facecolor='magenta')

ax.plot(dec7,major7,'k-')
ax.plot(dec7,minor7,'k-')
ax.fill_between(dec7,major7,minor7,facecolor='blue')


ax.set_xlim((-70., 30.0))
ax.set_ylim((0., 4.0))
ax.set_xlabel('Declination (degree)',size= 15)
ax.set_ylabel('beam (arcsec,100 GHz)', size = 15)


pl.show()

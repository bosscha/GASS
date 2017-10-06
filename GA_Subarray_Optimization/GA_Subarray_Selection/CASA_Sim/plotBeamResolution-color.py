#!/usr/bin/python

"""
display the minor/major beam resolution for different array configuration

"""

import sys
import pickle
import pylab as pl
import math

f = open('Resolution-C43-1.pickle', 'r')
data = pickle.load(f)
f.close()

dec0 = data[0]
minor0 = data[1]
major0 = data[2]
las0 = data[3]



f = open('Resolution-C43-2.pickle', 'r')
data = pickle.load(f)
f.close()

dec1 = data[0]
minor1 = data[1]
major1 = data[2]
las1 = data[3]


f = open('Resolution-C43-3.pickle', 'r')
data = pickle.load(f)
f.close()

dec2 = data[0]
minor2 = data[1]
major2 = data[2]
las2 = data[3]



f = open('Resolution-C43-4.pickle', 'r')
data = pickle.load(f)
f.close()

dec3 = data[0]
minor3 = data[1]
major3 = data[2]
las3 = data[3]



f = open('Resolution-C43-5.pickle', 'r')
data = pickle.load(f)
f.close()

dec4 = data[0]
minor4 = data[1]
major4 = data[2]
las4 = data[3]


f = open('Resolution-C43-6.pickle', 'r')
data = pickle.load(f)
f.close()

dec5 = data[0]
minor5 = data[1]
major5 = data[2]
las5 = data[3]


f = open('Resolution-C43-7.pickle', 'r')
data = pickle.load(f)
f.close()

dec6 = data[0]
minor6 = data[1]
major6 = data[2]
las6 = data[3]


f = open('Resolution-C43-8.pickle', 'r')
data = pickle.load(f)
f.close()

dec7 = data[0]
minor7 = data[1]
major7 = data[2]
las7 = data[3]


f = open('Resolution-C43-9.pickle', 'r')
data = pickle.load(f)
f.close()

dec8 = data[0]
minor8 = data[1]
major8 = data[2]
las8 = data[3]

f = open('Resolution-C43-10.pickle', 'r')
data = pickle.load(f)
f.close()

dec9 = data[0]
minor9 = data[1]
major9 = data[2]
las9 = data[3]

f = open('Resolution-ACA-std.pickle', 'r')
data = pickle.load(f)
f.close()

dec10 = data[0]
minor10 = data[1]
major10 = data[2]
las10 = data[3]


# print dec0[20]
# 
# print minor0[20]
# print major0[20]
# print minor1[20]
# print major1[20]
# 
# print minor2[20]
# print major2[20]
# print minor3[20]
# print major3[20]
# print minor4[20]
# print major4[20]
# 
# print minor5[20]
# print major5[20]
# 
# print minor6[20]
# print major6[20]
# print minor7[20]
# print major7[20]



# res0 = []
# res1 = []
# res2 = []
# res3 = []
# res4 = []
# res5 = []
# res6 = []
# res7 = []
# res8 = []
# 
# 
# 
# for i in range(0,len(minor0)):
#   res0.append(math.sqrt(minor0[i]*major0[i]))
#   res1.append(math.sqrt(minor1[i]*major1[i]))
#   res2.append(math.sqrt(minor2[i]*major2[i]))
#   res3.append(math.sqrt(minor3[i]*major3[i]))
#   res4.append(math.sqrt(minor4[i]*major4[i]))
#   res5.append(math.sqrt(minor5[i]*major5[i]))
#   res6.append(math.sqrt(minor6[i]*major6[i]))
#   res7.append(math.sqrt(minor7[i]*major7[i]))
#   res8.append(math.sqrt(minor8[i]*major8[i]))


fig = pl.figure()
ax = fig.add_subplot(111)

ax.semilogy(dec0,major0,'k-')
ax.semilogy(dec0,minor0,'k-')
ax.fill_between(dec0,major0,minor0,facecolor='yellow')

ax.semilogy(dec1,major1,'k-')
ax.semilogy(dec1,minor1,'k-')
ax.fill_between(dec1,major1,minor1,facecolor='blue')

ax.plot(dec2,major2,'k-')
ax.plot(dec2,minor2,'k-')
ax.fill_between(dec2,major2,minor2,facecolor='green')

ax.plot(dec3,major3,'k-')
ax.plot(dec3,minor3,'k-')
ax.fill_between(dec3,major3,minor3,facecolor='red')

ax.plot(dec3,major4,'k-')
ax.plot(dec3,minor4,'k-')
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


ax.plot(dec8,major8,'k-')
ax.plot(dec8,minor8,'k-')
ax.fill_between(dec8,major8,minor8,facecolor='green')


ax.plot(dec9,major9,'k-')
ax.plot(dec9,minor9,'k-')
ax.fill_between(dec9,major9,minor9,facecolor='red')

ax.plot(dec10,major10,'k-')
ax.plot(dec10,minor10,'k-')
ax.fill_between(dec10,major10,minor10,facecolor='yellow')



#ax.semilogy(dec0,res0,'b-')
#ax.semilogy(dec1,res1,'b-')
#ax.semilogy(dec2,res2,'b-')
#ax.semilogy(dec3,res3,'b-')
#ax.semilogy(dec4,res4,'b-')
#ax.semilogy(dec5,res5,'b-')
#ax.semilogy(dec6,res6,'b-')
#ax.semilogy(dec7,res7,'b-')
#ax.semilogy(dec8,res8,'r-')





ax.set_xlim((-70., 40.0))
ax.set_ylim((0.03, 30.0))

labels = [0.05, 0.1, 0.5, 1, 1.5, 2, 3, 4, 8, 20]
ax.set_yticks(labels)
ax.set_yticklabels([str(i) for i in labels])

ax.set_xlabel('Declination (degree)',size= 15)
ax.set_ylabel('Beam (arcsec,100 GHz)', size = 15)

ax.text(-20,19,'7-m',size = 15)
ax.text(-20,4.6,'C43',size = 15)

pl.savefig("beamResolutionDeclination-color-Cy5.png")

pl.show()

  






## plot the resolution for different Briggs factor
## C40-4.cfg // 100 GHz // 1h // -23d
## -2 Uniform
## 2 Natural


import numpy as np
import os
import math
import pickle
import pylab as pl


f = open("brigss.pickle","r")
data = pickle.load(f)
f.close()

briggs = data[0]
resolution = data[1]


fig = pl.figure()
ax = fig.add_subplot(111)

ax.plot(briggs,resolution,'b-')

c = 299792458.0
lam =  c / 100.0e9

lmax = 704.1 # this is for C40-4
rad2arcsec = 180.*3600. / 3.14159

k = rad2arcsec * lam / lmax

print k

xx = [-2.,2.]
yy = [k,k]

ax.plot(xx,yy,'k--')

ax.set_xlim((-2.,2.0))
ax.set_ylim((0.,1.6))
ax.text(-1.8,0.95,r'$\theta = \frac{\lambda}{L_{max}}$', size = 20)

ax.set_xlabel('Robust factor (Briggs)',size= 15)
ax.set_ylabel(r'$\theta$ (arcsec)', size = 15)

pl.savefig("briggs-C43-4.png")
pl.show()

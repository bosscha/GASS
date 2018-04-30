## To plot the UV radial distribution.

import UVW 
#s1 = 'sim/sim.ACA-9-02.ms'
s2 = 'cy5-1/sim/sim.ACA-std.ms' 
s3 = 'cy5-1/sim/sim.C43-1.ms'
s4 = 'cy5-1/cy5-1-full.int.ms' 

# k = 7.* 7. / (12.*12.)
#u1 = UVW.UVW(s1)
#r1, d1 = u1.radialDensity(0.,1000.,250)
#d1 = k * d1
#r1, c1 = u1.radialDensity(0.,1000.,250, norm = False)

u2 = UVW.UVW(s2)
r2, d2 = u2.radialDensity(0.,200.,150)
r2, c2 = u2.radialDensity(0.,200.,150, norm = False)
# u2.plotUVCoverage(-2000,2000,-2000,2000,'C36-3',saveFig=True,figfile='C36-3.uv.png')


u3 = UVW.UVW(s3)
r3, d3 = u3.radialDensity(0.,200.,150)
r3, c3 = u3.radialDensity(0.,200.,150, norm = False)
# u3.plotUVCoverage(-2000,2000,-2000,2000,'C36-6',saveFig=True,figfile='C36-6.uv.png')



u4 = UVW.UVW(s4)
r4, d4 = u4.radialDensity(0.,200.,150)
r4, c4 = u4.radialDensity(0.,200.,150, norm = False)
# u4.plotUVCoverage(-2000,2000,-2000,2000,'Full',saveFig=True,figfile='Full.uv.png')




fig = pl.figure()

ax = fig.add_subplot(111)
# ax.plot(r1, d1,'b-')
ax.plot(r2, d2,'k-')
ax.plot(r3, d3,'r-')
ax.plot(r4, d4,'g-')


# ax.plot(rad1,dens1[3],'c-')
ax.set_xlim((0., 200.0))
ax.set_ylim((0.,130))
# leg = ax.legend(('O-1-40','O-2-40','O-4-40','O-7-40','O-10-40','O-13-40','O-16-40','O-19-40'), 'upper right')
leg = ax.legend(('7-m','C43-1','Full'), 'upper right')
ax.grid(False)
ax.set_xlabel('UV radius (meter)')
ax.set_ylabel('weight')

pl.savefig("C43-1+ACA-uvdensity.png")
pl.show()

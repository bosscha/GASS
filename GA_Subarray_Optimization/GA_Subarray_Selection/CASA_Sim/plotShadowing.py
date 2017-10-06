import pylab as pl
import pickle


f = open('ACA-std-1h.pickle', 'r')
arr0 = pickle.load(f)
f.close()

f = open('C43-1-1h.pickle', 'r')
arr2 = pickle.load(f)
f.close()

f = open('C43-2-1h.pickle', 'r')
arr3 = pickle.load(f)
f.close()



fig = pl.figure()

ax = fig.add_subplot(111)

ax.plot(arr0[0],arr0[1],'ko-')
# ax.plot(arr1[0],arr1[1],'b*-')
ax.plot(arr2[0],arr2[1],'rd-')
ax.plot(arr3[0],arr3[1],'yD-')

#ax.plot(dec2,shad2,'b*-')
#ax.plot(dec3,shad3,'r+-')

ax.set_xlim((-90.,40.))
ax.set_ylim((0., 70.))



leg = ax.legend(('7-m','C43-1','C43-2',), 'upper middle')
ax.grid(False)
ax.set_xlabel('Declination (Degree)')
ax.set_ylabel('Shadowing (%)')

pl.savefig("shadowing-Cy5.png")
pl.show()
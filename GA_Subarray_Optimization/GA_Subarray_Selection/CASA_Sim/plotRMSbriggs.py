import pickle
import pylab as pl



f = open('c40-1.rms.brigss.pickle')
datac401 = pickle.load(f)
f.close()

f = open('c40-2.rms.brigss.pickle')
datac402 = pickle.load(f)
f.close()


f = open('c40-3.rms.brigss.pickle')
datac403 = pickle.load(f)
f.close()

f = open('c40-9.rms.brigss.pickle')
datac409 = pickle.load(f)
f.close()


x1 = datac401[0]
y1 = datac401[1]
print("Factor C40-1: %f"%(y1[62]/y1[-1]))

x2 = datac402[0]
y2 = datac402[1]
print("Factor C40-2: %f"%(y2[62]/y2[-1]))
      
x3 = datac403[0]
y3 = datac403[1]
print("Factor C40-3: %f"%(y3[62]/y3[-1]))

x9 = datac409[0]
y9 = datac409[1]
print("Factor C40-9: %f"%(y9[62]/y9[-1]))


# print("Factor RMS : %f \n"%(yy[62]/yy[-1]))

fig = pl.figure()
ax = fig.add_subplot(111)

ax.plot(x1, y1 /y1[-1],'k-')
ax.plot(x2, y2 /y2[-1],'g-')
ax.plot(x3, y3 /y3[-1],'r-')
ax.plot(x9, y9 /y9[-1],'c-')


#ax.set_xlim((-90.,40.))
#ax.set_ylim((0., 70.))



leg = ax.legend(('C40-1','C40-2','C40-3','C40-9'), 'upper right')


xx0 =[-2,2]
yy0 = [1.,1.]
ax.plot(xx0,yy0,'y--')

pl.show()




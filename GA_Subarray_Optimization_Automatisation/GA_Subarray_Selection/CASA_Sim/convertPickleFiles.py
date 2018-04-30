import glob
import math
import pickle

fnames = glob.glob('*pickle')

for i in range(len(fnames)):

    f = open(fnames[i])
    fc = pickle.load(f)
    f.close()

    f = open(fnames[i]+'.dat', 'w')
    
    print >> f, '# Declination  Resolution  MRS '
    print >> f, '#  (degree)      (arcsec)  (arcsec) '

    for j in range(len(fc[0])):
        print >> f, "%4d   %.3f   %.3f  " %(fc[0][j], math.sqrt(fc[1][j]*fc[2][j]), fc[3][j])

    f.close()

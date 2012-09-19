import numpy as np
import scipy as sp
import pylab as pl
import scipy.optimize as opt

def MB(v,T):
    return 2*v/T * np.exp(-v**2/T)

def fit(T, vreal, preal):
    pguess = MB(vreal,T)
    return ((pguess - preal)**2).sum()

def showfit(T, vreal, preal):
    pl.figure()
    pl.plot(vreal, preal, 'o', label='Simulation')
    pl.plot(vreal, MB(vreal,T), label='2D MB')
    pl.legend()
    pl.xlabel(r'$|r|$', fontsize=20)
    pl.ylabel(r'$P(r)$', fontsize=20)
    pl.title("Velocity Distribution in Moshpit", fontsize=20)
    pl.show()

if __name__ == "__main__":
    r = np.fromfile(open("velocities.txt", "rb"))
    r = r[r<6]
    
    h = pl.hist(r, bins=80)
    pl.show()
    
    vreal = h[1][:-1]
    preal = 1.*h[0] / h[0].sum() / (vreal[1] - vreal[0]) 
    
    f = opt.fmin(fit, [6], args=(vreal,preal), xtol=1e-8)
    showfit(f, vreal, preal)

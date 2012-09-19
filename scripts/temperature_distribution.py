from numpy import *
from scipy import *
from pylab import *
import scipy.optimize as opt

vreal = arange(0, 3, 3.0/50)

def dofit(damp, rad):
    r = loadtxt("temp_%0.2f" % damp + ".txt")
    p = r[rad,:]
    preal = 1.*p / p.sum() / (vreal[1] - vreal[0]) 
    f = opt.fmin(fit, [3], args=(vreal,preal), xtol=1e-8)
    return f
    #showfit(f, vreal,preal)

def fit(T, vreal, preal):
    pguess = 2*vreal/T * exp(-vreal**2/T)
    return ((pguess - preal)**2).sum()

def showfit(T, vreal, preal):
    plot(vreal, preal, 'o', label='Simulation')
    plot(vreal, 2*vreal/T * exp(-vreal**2/T), label='2D MB')

figure()
for i in arange(0.05, 0.51, 0.05):
    temps = []
    for j in range(10):
        temps.append(dofit(i,j))
    plot(range(len(temps)), temps, 'o-', label=r"$\beta=%0.2f$" % i)

xlabel(r'$|r|$', fontsize=20)
ylabel(r'$T(r)$', fontsize=20)
title("Temperature Distribution in Moshpit", fontsize=20)
#legend()
savefig("temperature.png")

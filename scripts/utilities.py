import numpy as np
import scipy as sp
import pylab as pl
import scipy.optimize as opt
import re, shutil, os

from subprocess import Popen, PIPE, STDOUT
import time

#===============================================
# utilities to run the simulation
#===============================================
def changeLine(match, replacement, filename):
    newfl = filename+"~"
    with open(newfl, "w") as newfile:
        contents = open(filename).readlines()
        for line in contents:
            if len(re.findall(match, line)) > 0:
                line = replacement+"\n"
            newfile.write(line)
    shutil.move(newfl, filename)

def setOptions(fps=0, opengl=0, velocities=0, temperature=0, timeseries=0):
    cvel = "" if velocities else "//"
    ctmp = "" if temperature else "//"
    cang = "" if timeseries else "//"
    dop  = "1" if opengl else "0"
    dof  = "1" if fps else "0"
    changeLine(r"#define VELOCITY_DISTRIBUTION", cvel+"#define VELOCITY_DISTRIBUTION", "../main.c")
    changeLine(r"#define TEMPERATURE_BINS",      ctmp+"#define TEMPERATURE_BINS",      "../main.c")
    changeLine(r"#define ANGULARMOM_TIMESERIES", cang+"#define ANGULARMOM_TIMESERIES", "../main.c")
    changeLine(r"^DOPLOT", "DOPLOT = "+dop, "../Makefile")
    changeLine(r"^FPS",    "FPS    = "+dof, "../Makefile")
    print "Now building entbody"
    os.system("cd .. && make clean && make && cd -")

def launchSingleMoshpit(alpha, eta, seed, damp=1.0):
    return Popen("nice -n 20 ../entbody "+str(alpha)+" "+str(eta)+" "+str(seed)+" "+str(damp), 
            shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)

def runSingleMoshpit(alpha, eta, seed, damp=1.0):
    proc = [launchSingleMoshpit(alpha, eta, seed, damp)]
    while proc:
        for p in proc:
            retcode = p.poll()
            if retcode is not None:
                out = p.stdout.readlines()[-1]
                proc.remove(p)
    return [float(o) for o in out.split(' ')]

def runEntireSlice(samples=50, clump=10, alpha_range=(0.0,4.0,4.0/60), filename="runs.txt"):
    setOptions()
    curr = 0
    file = open(filename, "w", 0)
    
    for alpha in np.arange(*alpha_range):
        for eta in np.arange(0.0, 3.0, 0.06):
            start = time.time()
            data = [[] for i in xrange(12)]
            for a in range(0, samples, clump):
                curr += clump 
                procs = [launchSingleMoshpit(alpha,eta,seed) for seed in range(curr, curr+clump)]
                while procs:
                    for p in procs:
                        retcode = p.poll()
                        if retcode is not None:
                            out = p.stdout.read().split(' ')
                            for i in xrange(len(out)):
                                data[i].append(float(out[i]))
                            procs.remove(p)
            strout = str(alpha)+" "+str(eta)+" "
            strout += " ".join([str(np.mean(d)) for d in data])+" "
            strout += " ".join([str(np.std(d)) for d in data])
            strout += "\n"
            file.write(strout)
            end = time.time()
            print alpha, " ", eta, " ", end - start

#=====================================================
# helper functions that generate a MB fit
#=====================================================
def MB(v,T):
    return 2*v/T * np.exp(-v**2/T)

def fitToMB(T, vreal, preal):
    pguess = MB(vreal,T)
    return ((pguess - preal)**2).sum()

def showFitToMB(T, vreal, preal):
    pl.figure()
    pl.plot(vreal, preal, 'o', label='Simulation')
    pl.plot(vreal, MB(vreal,T), label='2D MB')
    pl.legend()
    pl.xlabel(r'$|r|$', fontsize=20)
    pl.ylabel(r'$P(r)$', fontsize=20)
    pl.title("Velocity Distribution in Moshpit", fontsize=20)
    pl.savefig("velocitydist.png")

def runVelocityFit():
    setOptions(velocities=1)
    runSingleMoshpit(0.2,0.6,1,0.3)
    setOptions(velocities=1, opengl=1)

    r = np.fromfile(open("velocities.txt", "rb"))
    r = r[r<6]
    
    h = pl.hist(r, bins=80)
    pl.show()
    
    vreal = h[1][:-1]
    preal = 1.*h[0] / h[0].sum() / (vreal[1] - vreal[0]) 
    
    f = opt.fmin(fitToMB, [6], args=(vreal,preal), xtol=1e-8, disp=0)
    showFitToMB(f, vreal, preal)


#====================================================
# functions that generate the temperature profile
#====================================================
def runTemperatureFits():
    def dofit(damp, rad):
        r = np.loadtxt("temp_%0.2f" % damp + ".txt")
        p = r[rad,:]
        vreal = np.arange(0, 125, 125.0/50)
        preal = 1.*p / p.sum() / (vreal[1] - vreal[0]) 
        f = opt.fmin(fitToMB, [30], args=(vreal,preal), xtol=1e-8, disp=0)
        return f[0]

    setOptions(temperature=1)
    pl.figure()
    for i in np.arange(0.1, 0.5, 0.1):
        runSingleMoshpit(0.0,0.0,1,i)
        temps = []
        for j in range(6):
            temps.append(dofit(i,j))
        pl.plot(range(len(temps)), temps, 'o-', label=r"$\beta=%0.2f$" % i)
    setOptions(opengl=1)

    pl.xlabel(r'$|r|$', fontsize=20)
    pl.ylabel(r'$T(r)$', fontsize=20)
    pl.title("Temperature Distribution in Moshpit", fontsize=20)
    pl.savefig("temperature.png")

def runTemperatureSlice(beta=0.25):
    def dofit(damp, rad):
        r = np.loadtxt("temp_%0.2f" % damp + ".txt")
        p = r[rad,:]
        vreal = np.arange(0, 125, 125.0/50)
        preal = 1.*p / p.sum() / (vreal[1] - vreal[0]) 
        f = opt.fmin(fitToMB, [30], args=(vreal,preal), xtol=1e-8, disp=0)
        pl.plot(preal, 'o', label=str(rad))
        pl.plot(MB(vreal, f[0]), '-')
        return f[0]

    setOptions(temperature=1)
    pl.figure()
    runSingleMoshpit(0.1,1.9,1,beta)
    temps = []
    for j in range(10):
        ttemp = dofit(beta,j)
        #pl.plot(range(len(temps)), temps, 'o-', label=r"$\beta=%0.2f$" % i)
    pl.legend()
    setOptions(opengl=1)

    pl.xlabel(r'$|r|$', fontsize=20)
    pl.ylabel(r'$T(r)$', fontsize=20)
    pl.title("Temperature Distribution in Moshpit", fontsize=20)
    pl.savefig("temperature.png")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "all":
        print "This is going to take a while..."
        runVelocityFit()
        runTemperatureFits()
        runEntireSlice()

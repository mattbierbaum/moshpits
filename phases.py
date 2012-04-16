from numpy import *
from scipy import *
from pylab import *

from subprocess import Popen, PIPE, STDOUT
import time

samples = 200
curr = 0
vorticies = []
file = open("runs.txt", "w", 0)

for alpha in arange(0.0, 1.0, 0.02):
    for eta in arange(0.0, 10.0, 0.2):
        start = time.time()
        vorts = []
        curr += samples
        procs = [Popen("nice -n 20 ./entbody "+str(alpha)+" "+str(eta)+" "+str(seed), shell=True, stdin=PIPE, stdout=PIPE, close_fds=True) for seed in range(curr, curr+samples)]
        while procs:
            for p in procs:
                retcode = p.poll()
                if retcode is not None:
                    vorts.append(abs(float(p.stdout.read())))
                    procs.remove(p)
        vorts = array(vorts)
        vorticies.append([alpha, eta,vorts.mean(),vorts.std()])
        file.write(str(alpha)+" "+str(eta)+" "+str(vorts.mean())+" "+str(vorts.std())+"\n")
        end = time.time()
        print alpha, " ", eta, " ", end - start

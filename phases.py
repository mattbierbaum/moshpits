from numpy import *
from scipy import *
from pylab import *

from subprocess import Popen, PIPE, STDOUT
import time

samples = 50 
curr = 0
vorticies = []
file = open("runs.txt", "w", 0)

for alpha in arange(0.0, 5.0, 0.1):
    for eta in arange(0.0, 5.0, 0.1):
        start = time.time()
        avg,avgs,sq,sqs = [], [], [], []
        curr += samples
        procs = [Popen("nice -n 20 ./entbody "+str(alpha)+" "+str(eta)+" "+str(seed), shell=True, stdin=PIPE, stdout=PIPE, close_fds=True) for seed in range(curr, curr+samples)]
        while procs:
            for p in procs:
                retcode = p.poll()
                if retcode is not None:
                    tavg, tavgs, tsq, tsqs = [float(s) for s in  p.stdout.read().split(' ')]
                    avg.append(tavg)
                    avgs.append(tavgs)
                    sq.append(tsq)
                    sqs.append(tsqs)
                    procs.remove(p)
        avg  = array(avg)
        avgs = array(avgs)
        sq   = array(sq)
        sqs  = array(sqs) 
        file.write(str(alpha)+" "+str(eta)+" "+str(avg.mean())+" "+str(avgs.mean())+" "+str(sq.mean())+" "+str(sqs.mean())+" "+str(avg.std())+" "+str(avgs.std())+" "+str(sq.std())+" "+str(sqs.std())+"\n")
        end = time.time()
        print alpha, " ", eta, " ", end - start

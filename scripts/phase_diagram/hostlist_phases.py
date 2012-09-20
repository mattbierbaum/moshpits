from numpy import *
from scipy import *
from pylab import *

from subprocess import Popen, PIPE, STDOUT
import time

from socket import gethostname
host = gethostname()
hostn = int(host[3]+host[4])
sa = 4.0/20. * (hostn-1)
se = 4.0/20. * (hostn)
ss = 4.0/60.

samples = 200 
curr = 0
vorticies = []
file = open("runs.txt.%s" % host, "w", 0)

for alpha in arange(sa, se, ss):
    for eta in arange(0.0, 3.0, 0.06):
        start = time.time()
        avg,avgs,sq,sqs = [], [], [], []
        mxa,mxs,mya,mys = [], [], [], []
        mxsa, mxss, mysa, myss = [], [], [], []
        curr += samples
        procs = [Popen("nice -n 20 ./entbody "+str(alpha)+" "+str(eta)+" "+str(seed), shell=True, stdin=PIPE, stdout=PIPE, close_fds=True) for seed in range(curr, curr+samples)]
        while procs:
            for p in procs:
                retcode = p.poll()
                if retcode is not None:
                    tavg, tavgs, tsq, tsqs, tmxa,tmxs,tmya,tmys, tmxsa, tmxss, tmysa, tmyss = [float(s) for s in  p.stdout.read().split(' ')]
                    avg.append(tavg)
                    avgs.append(tavgs)
                    sq.append(tsq)
                    sqs.append(tsqs)
                    mxa.append(tmxa)
                    mxs.append(tmxs)
                    mya.append(tmya)
                    mys.append(tmys)
                    mxsa.append(tmxsa)
                    mxss.append(tmxss)
                    mysa.append(tmysa)
                    myss.append(tmyss)
                    procs.remove(p)
        avg  = array(avg)
        avgs = array(avgs)
        sq   = array(sq)
        sqs  = array(sqs)
        mxa  = array(mxa)
        mxs  = array(mxs)
        mya  = array(mya)
        mys  = array(mys)
        mxsa = array(mxsa)
        mxss = array(mxss)
        mysa = array(mysa)
        myss = array(myss)
        file.write(str(alpha)+" "+str(eta)+" "+str(avg.mean())+" "+str(avgs.mean())+" "+str(sq.mean())+" "+str(sqs.mean())+" "+str(avg.std())+" "+str(avgs.std())+" "+str(sq.std())+" "+str(sqs.std())+" "+str(mxa.mean())+" "+str(mya.mean())+" "+str(mxs.mean())+" "+str(mys.mean())+" "+str(mxa.std())+" "+str(mya.std())+" "+str(mxsa.mean())+" "+str(mxss.mean())+" "+str(mysa.mean())+" "+str(myss.mean())+" "+str(mxsa.std())+" "+str(mysa.std())+"\n")
        end = time.time()
        print alpha, " ", eta, " ", end - start

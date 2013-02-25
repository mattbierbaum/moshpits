from utilities import runEntireSlice
from socket import gethostname

host = gethostname()
hostn = int(host[3]+host[4])
sa = 4.0/20. * (hostn-1)
se = 4.0/20. * (hostn)
ss = 4.0/60.

samples = 200 
clump   = 10

runEntireSlice(samples, clump, (sa,se,ss), "runs.txt.%s"%host)

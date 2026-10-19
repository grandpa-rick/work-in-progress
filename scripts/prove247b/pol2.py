import scert,sys
from scert import *
scert.LIMS=(0,)
N=int(sys.argv[1]); bad=[]
for lam,mu,nu,a,b in allA(N):
    if not start('A',N,lam,mu,nu,a,b): bad.append((lam,mu,nu,a,b))
print(len(bad)); 
for x in bad: print(x)

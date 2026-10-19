from scert import *
import sys,time
N=int(sys.argv[1])
t=time.time(); bad=[]; tot=0
for lam,mu,nu,a,b in allA(N):
    tot+=1
    if not start('A',N,lam,mu,nu,a,b): bad.append((lam,mu,nu,a,b))
print(N,'allowed',tot,'uncert',len(bad),'time',round(time.time()-t),STATS); print(bad[:10])

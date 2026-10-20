import sys,time,collections; sys.path.insert(0,"../prove247b"); sys.path.insert(0,'../prove247')
from scert5 import *
N=int(sys.argv[1]); t=time.time(); bad=[]; tot=0
for lam,mu,nu,a,b in allA(N):
    tot+=1
    if not start('A',N,lam,mu,nu,a,b): bad.append((lam,mu,nu,a,b))
print('N',N,'allowed',tot,'uncert',len(bad),bad[:10],round(time.time()-t),'s',PSTAT,flush=True)
print(collections.Counter(x[2] for x in bad).most_common(10))

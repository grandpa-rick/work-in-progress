import sys,ast; sys.path.insert(0,'../prove247')
import mp
from scert import *
N=7; bad=[]
for lam,mu,nu,a,b in allA(N):
    if not start('A',N,lam,mu,nu,a,b): bad.append((lam,mu,nu,a,b))
print('uncert',len(bad),flush=True)
sy={}
z=0
for lam,mu,nu,a,b in bad:
    if (a,b) not in sy: sy[a,b]=mp.Sys(N,5,a,b)
    for c in (11,12,13):
        if mp.A(sy[a,b],a,b,lam,mu,nu,c)==0: z+=1; print('ZERO?',lam,mu,nu,a,b,c,flush=True)
print('checked',len(bad),'x3 c; zeros',z)
import collections; print(collections.Counter((x[3],x[4]) for x in bad)); print(collections.Counter(x[2] for x in bad).most_common(8))

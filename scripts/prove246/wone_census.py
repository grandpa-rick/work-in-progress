import sys, itertools, time; sys.path.insert(0,'.')
import numpy as np
from wone import *
N=int(sys.argv[1]); wt=int(sys.argv[2]); pairs=[tuple(map(int,p.split(','))) for p in sys.argv[3:]] or list(itertools.combinations(range(1,N+1),2))
Ss=subsets(N)
tri=[(l,m,nu) for i,l in enumerate(Ss) for m in Ss[i:] for nu in Ss if geq(nu,l) and geq(nu,m) and len(nu)-max(len(l),len(m))>=2]
def allowed(l,m,nu,a,b):
    rest=[c for c in range(1,N+1) if c not in (a,b)]
    for k in range(len(rest)+1):
        for Mpp in itertools.combinations(rest,k):
            M=tuple(sorted(Mpp+(a,b)))
            if geq(nu,M) and geq(Mpp,l) and geq(Mpp,m): return True
    return False
cs=list(range(-12,13,2)); t0=time.time(); summ={}
for a,b in pairs:
    Bs={c:None for c in cs}
    B=BW(N,a,b,wt,Fr(3),Fr(1,10**10)); B2=BW(N,a,b,wt,Fr(3),Fr(1,10**7))
    for l,m,nu in tri:
        if not allowed(l,m,nu,a,b): continue
        v=np.array([float(B.residc(l,m,nu,c)) for c in cs]); v2=np.array([float(B2.residc(l,m,nu,c)) for c in cs])
        conv=np.max(np.abs(v-v2))/max(1e-30,np.max(np.abs(v)))
        X=np.array(cs,dtype=float)
        for d in range(0,8):
            co=np.polyfit(X,v,d); err=np.max(np.abs(np.polyval(co,X)-v))/max(1e-30,np.max(np.abs(v)))
            if err<1e-6: break
        roots=np.roots(co) if d>0 else []
        rr=sorted([round(z.real,3) if abs(z.imag)<1e-4 else complex(round(z.real,3),round(z.imag,3)) for z in roots],key=lambda z:(abs(z),str(z)))
        key=(d,str(rr)); summ[key]=summ.get(key,0)+1
        print((a,b),l,m,nu,'deg',d,'roots',rr,'conv %.1e'%conv, 'val(-2N)=%.3e'%np.polyval(co,-2*N),flush=True)
print('SUMMARY',summ,'%.0fs'%(time.time()-t0))

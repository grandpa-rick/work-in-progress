"""Clio UID 419 ask 6: do the a,b-not-in-nu cases of E != empty give actual H_ab residues/poles?
Same engine and detector as ../prove246/census.py (resid at seeds 3,11; wt0 = q=0 basis, wt1 = q=inf basis).
Reports, per N, the allowed cases with {a,b} not a subset of nu, and how many have R_inf != 0 (and R_0 != 0)."""
import sys, itertools; sys.path.insert(0,'../prove246')
from fast import *
N=int(sys.argv[1])
Ss=subsets(N)
tri=[(l,m,nu) for i,l in enumerate(Ss) for m in Ss[i:] for nu in Ss if geq(nu,l) and geq(nu,m) and len(nu)-max(len(l),len(m))>=2]
def E(l,m,nu,a,b):
    rest=[c for c in range(1,N+1) if c not in (a,b)]
    return [Mpp for k in range(len(rest)+1) for Mpp in itertools.combinations(rest,k)
            if geq(nu,tuple(sorted(Mpp+(a,b)))) and geq(Mpp,l) and geq(Mpp,m)]
tot=0; nzi=0; nz0=0; ex=[]
for a,b in itertools.combinations(range(1,N+1),2):
    B={(wt,sd):Basis(N,a,b,sd,wt) for wt in (0,1) for sd in (3,11)}
    for l,m,nu in tri:
        if a in nu and b in nu: continue
        e=E(l,m,nu,a,b)
        if not e: continue
        tot+=1
        zi=all(B[1,sd].resid(l,m,nu)==0 for sd in (3,11))
        z0=all(B[0,sd].resid(l,m,nu)==0 for sd in (3,11))
        nzi+= (not zi); nz0+=(not z0)
        if (l,m,nu,a,b)==((),(1,),(1,2,4),1,3) or len(ex)<3: ex.append(((l,m,nu,a,b),e,'Rinf!=0' if not zi else 'Rinf=0','R0!=0' if not z0 else 'R0=0'))
print('N',N,'allowed cases with {a,b} not in nu:',tot,' R_inf!=0:',nzi,' R_0!=0:',nz0)
for x in ex: print('  ',x)

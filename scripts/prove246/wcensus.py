"""R(c) = Res l^F_nu(1_{>=ab} * prod_{m in M} w_m^c * F_lam F_mu), F = phi0 (wt0) or N-free phi^inf form (wt1: exponent 2r)."""
import sys, itertools, time; sys.path.insert(0,'.')
import fast
from fast import *
N=int(sys.argv[1]); cs=list(map(int,sys.argv[2].split(','))); pairs=[tuple(map(int,p.split(','))) for p in sys.argv[3:]] or list(itertools.combinations(range(1,N+1),2))
Ss=subsets(N)
tri=[(l,m,nu) for i,l in enumerate(Ss) for m in Ss[i:] for nu in Ss if geq(nu,l) and geq(nu,m) and len(nu)-max(len(l),len(m))>=2]
def allowed(l,m,nu,a,b):
    rest=[c for c in range(1,N+1) if c not in (a,b)]
    for k in range(len(rest)+1):
        for Mpp in itertools.combinations(rest,k):
            M=tuple(sorted(Mpp+(a,b)))
            if geq(nu,M) and geq(Mpp,l) and geq(Mpp,m): return True
    return False
class B2(Basis):
    def A(s,M,I):
        if (M,I) not in s.cache: s.cache[M,I]=Fv(I,M,s.w,s.wt,NN=0) if geq(M,I) else S.const(0)
        return s.cache[M,I]
    def residc(s,lam,mu,nu,c):
        a,b=s.a,s.b; f={}
        for M in subsets(s.N):
            if a in M and b in M and geq(nu,M):
                wgt=S.const(1)
                for m in M: wgt=wgt*(s.w[m-1]**c)
                f[M]=s.A(M,lam)*s.A(M,mu)*wgt
        return s.coords(f,nu).coef(-1)
t0=time.time()
for wt in (0,1):
  for c in cs:
    fails=[]; tot=0
    for a,b in pairs:
        Bs=[B2(N,a,b,sd,wt) for sd in (3,11)]
        for l,m,nu in tri:
            if not allowed(l,m,nu,a,b): continue
            tot+=1
            if all(B.residc(l,m,nu,c)==0 for B in Bs): fails.append((a,b,l,m,nu))
    print('wt',wt,'c',c,'allowed',tot,'zeros',len(fails),fails[:6],'%.0fs'%(time.time()-t0),flush=True)

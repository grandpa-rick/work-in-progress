# w->1 limit of R^A on H_ab: polynomial in c
import mpq, itertools, sys, numpy as np
from mpq import *
from fractions import Fraction as Fr
class SysQ(Sys):
    def __init__(s,n,dl,al,a,b):
        s.n=n; s.w=[S.const(1+al[i]*dl) for i in range(n)]
        s.w[b-1]=S(0,[1/(1+al[a-1]*dl)]*K)   # w_b = 1/(w_a (1-eps))
        s.cache={}; s.subs=subsets(n); s.down={}
def RA(N,a,b,lam,mu,nu,c,al,dl):
    s=SysQ(N,dl,al,a,b); return A(s,a,b,lam,mu,nu,c)
def lead(N,a,b,lam,mu,nu,c,al,norm,d=None):
    vs=[]
    for k in (1,2,4):
        dl=Fr(k,10**14); vs.append(RA(N,a,b,lam,mu,nu,c,al,dl)/RA(N,*norm,al=al,dl=dl))
    import math
    if d is None:
        if vs[0]==0: return None,0.0
        d=round(math.log2(abs(float(vs[1]/vs[0]))))
    g=[v/(Fr(k,10**14))**d for v,k in zip(vs,(1,2,4))]
    return d,float((8*g[0]-6*g[1]+g[2])/3)
if __name__=="__main__":
    N=int(sys.argv[1]); import random; random.seed(int(sys.argv[2]))
    al=[Fr(random.randint(1,20),random.randint(1,20))*random.choice([1,-1]) for _ in range(N)]
    Uall=list(range(1,N+1)); subs=[tuple(S_) for r in range(N+1) for S_ in itertools.combinations(Uall,r)]
    dim=N*(N+1)//2; Umax=2*N-4
    norm=(N-1,N,(),(),tuple(Uall),2*N-2)
    cs=list(range(Umax+1-2*(dim+2),Umax+2+2*(dim+2)))
    pairs=eval(sys.argv[3]) if len(sys.argv)>3 else [(a,b) for a in Uall for b in Uall if a<b]
    for (a,b) in pairs:
     for nu in subs:
      for lam in subs:
       for mu in subs:
        if lam>mu: continue
        if not EA(N,a,b,lam,mu,nu): continue
        d,_=lead(N,a,b,lam,mu,nu,Umax+7,al,norm)
        if d is None: print('ZERO at c=U+7??',a,b,lam,mu,nu); continue
        vals=[lead(N,a,b,lam,mu,nu,c,al,norm,d)[1] for c in cs]
        co=np.polyfit(cs,vals,dim+2); sc=max(abs(v) for v in vals)
        co=np.trim_zeros(np.where(abs(co)*max(abs(np.array(cs)))**np.arange(len(co))[::-1]<1e-7*sc,0,co),'f')
        r=np.roots(co) if len(co)>1 else []
        mr=max([x.real for x in r],default=None)
        intabove=[c for c in range(Umax+1,Umax+40) if abs(np.polyval(co,c))<1e-6*sc]
        print((a,b),lam,mu,nu,'d',d,'deg',len(co)-1,'maxRe',None if mr is None else round(mr,3),'intzeros>U',intabove,flush=True)

import sys, itertools; sys.path.insert(0,'.')
import sympy as sp
from fast import subsets, geq
from cf import T0
def Fs(kappa,M,w,wt,N):
    r=len(kappa); tot=0
    if not geq(M,kappa): return sp.Integer(0)
    for h in itertools.permutations(M,r):
        if any(h[j]<kappa[j] for j in range(r)): continue
        t=T0(kappa,h,M,w)
        if wt: t*=sp.prod([w[h[j]-1]**(kappa[j]-2*(j+1)) for j in range(r)])
        tot+=t
    if wt: tot*=sp.prod([w[m-1]**(2*r-2*N) for m in M])
    return sp.factor(sp.cancel(tot))
def main():
  pass
if __name__=="__main__":
 N=int(sys.argv[1]) ; a,b=int(sys.argv[2]),int(sys.argv[3]); nu=tuple(map(int,sys.argv[4]))
 wt=int(sys.argv[5])
 w=list(sp.symbols('w1:%d'%(N+1))); x,eps=sp.symbols('x eps')
 down=sorted([M for M in subsets(N) if geq(nu,M)],key=lambda T:(len(T),sum(T)))
 sub={w[b-1]:1/(x*(1-eps)), w[a-1]:x}
 A={(M,I):sp.cancel(Fs(I,M,w,wt,N).subs(sub)) for M in down for I in down}
 def coords(f):
     c={}
     for M in down:
         v=f.get(M,0)-sum(c[I]*A[M,I] for I in c if A[M,I]!=0)
         c[M]=sp.cancel(v/A[M,M])
     return c
 for M in down:
     if a in M and b in M:
         r=sp.factor(sp.cancel(sp.cancel(coords({M:1})[nu]*eps).subs(eps,0)))
         print(M, r, flush=True)

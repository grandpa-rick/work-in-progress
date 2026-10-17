import sys, sympy as sp, time
from engine import q, subsets, geq, W, is_laurent
X=sys.argv[1]; N=int(sys.argv[2])
S=subsets(N); t0=time.time()
A={}
for M in S:
  for I in S:
    A[M,I] = W(X,N,I,M,sp.Integer(0))[0] if geq(M,I) else sp.Integer(0)
w=sp.symbols('w1:%d'%(N+1))
print('A0 built',time.time()-t0,flush=True)
order=sorted(S,key=lambda T:(len(T),sum(T)))
bad=0; nz2=[]
for l in S:
  for m in S:
    if order.index(m)<order.index(l): continue
    c={}
    for nu in order:
      rhs=A[nu,l]*A[nu,m]-sum(c[k]*A[nu,k] for k in c if A[nu,k]!=0)
      c[nu]=sp.factor(sp.cancel(rhs/A[nu,nu]))
    for nu,v in c.items():
      if v==0: continue
      if not is_laurent(v,w): bad+=1; print('NONLAURENT',l,m,nu,v)
      jump=len(nu)-max(len(l),len(m))
      if jump>=2: nz2.append((l,m,nu,v)); print('JUMP>=2 nonzero c0:',l,m,nu,v,flush=True)
print('nonlaurent',bad,'jump>=2 nonzero',len(nz2),'time',time.time()-t0)

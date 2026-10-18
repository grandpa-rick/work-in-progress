"""Explore q=infinity leading form of Ã vs w^{-1}-image of q=0 form. N from argv."""
import sys, pickle, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import q, subsets, geq
N=int(sys.argv[1]); S=subsets(N); w=sp.symbols('w1:%d'%(N+1)); dim=N*(N+1)//2
A=pickle.load(open('/home/agent/projects/work-in-progress/scripts/wake245/A_full_N%d.pkl'%N,'rb'))
p=sp.Symbol('p'); inv={x:1/x for x in w}
for I in S:
  for M in S:
    a=A[M,I]
    if a==0: continue
    n,d=sp.fraction(sp.cancel(a)); deg=sp.Poly(n,q).degree()-sp.Poly(d,q).degree()
    finf=sp.factor(sp.cancel(sp.cancel(a.subs(q,1/p)*p**deg).subs(p,0)))
    f0=sp.factor(a.subs(q,0))
    r=sp.factor(sp.cancel(finf/f0.subs(inv,simultaneous=True)))
    print('I=%s M=%s deg=%d (dim-codimI=%d) phi_inf=%s | ratio phi_inf/phi0(w^-1)=%s'%(I,M,deg,dim-sum(I),finf,r))

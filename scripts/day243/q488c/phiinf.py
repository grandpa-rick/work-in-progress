import sympy as sp
from engine import subsets, geq
from fullc import A, q
N=3; w=list(sp.symbols('w1:4')); t=sp.Symbol('t'); p=sp.Symbol('p'); d=N*(N+1)//2
e={1:1,2:3,3:9}
for M in subsets(N):
  row=[]
  for I in subsets(N):
    if not geq(M,I): row.append('.'); continue
    a=A(N,M,I,w); k=d-sum(I)
    fi=sp.cancel(a.subs(q,1/p)*p**k).subs(p,0)
    ft=sp.cancel(fi.subs({w[i-1]:t**(-e[i]) for i in e}))
    n_,dd=sp.fraction(ft); Pn=sp.Poly(n_,t); Pd=sp.Poly(dd,t)
    ln=min(Pn.monoms())[0]; ld=min(Pd.monoms())[0]
    cn=Pn.coeff_monomial(t**ln); cd=Pd.coeff_monomial(t**ld)
    row.append('%s*t^%d'%(sp.nsimplify(cn/cd),ln-ld))
  print(M,row)

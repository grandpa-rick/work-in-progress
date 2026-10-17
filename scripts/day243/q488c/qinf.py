import sympy as sp, itertools
from fullc import A, q
from engine import subsets, geq
N=2; w=list(sp.symbols('w1:3'))
for M in subsets(N):
  for I in subsets(N):
    if not geq(M,I): continue
    a=A(N,M,I,w)
    p=sp.Symbol('p')
    ai=sp.factor(sp.cancel(a.subs(q,1/p)))
    print(M,I,'A(q)=',sp.factor(a),'   A(1/q)=',ai)

import sys, sympy as sp
from engine import q, subsets, geq, W
X=sys.argv[1]; N=int(sys.argv[2])
S=subsets(N)
for M in S:
  for I in S:
    if not geq(M,I): continue
    a,w=W(X,N,I,M,q)
    a0=sp.factor(a.subs(q,0)); a1=sp.factor(sp.diff(a,q).subs(q,0))
    print(M,I,' A0=',a0,'   A1=',a1)

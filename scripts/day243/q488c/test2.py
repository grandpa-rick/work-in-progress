import sympy as sp
from engine import subsets, geq, W, q
from engine2 import A01
for N in (2,3):
  for M in subsets(N):
    for I in subsets(N):
      if not geq(M,I): continue
      a,w=W('C',N,I,M,q)
      b0,b1=A01(N,M,I)
      ok = sp.cancel(a.subs(q,0)-b0)==0 and sp.cancel(sp.diff(a,q).subs(q,0)-b1)==0
      if not ok: print('MISMATCH',N,M,I)
  print('N',N,'done')

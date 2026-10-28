import sys, sympy as sp
sys.path.insert(0,'prove262'); from q46 import *
m=3; k=(3,1); R=[r for r in strict(m) if sub(r,k)]; a=[3,2,1]
M=sp.Matrix([[sp.expand(spec(chw(r,j,m),a)) for r in R] for j in range(8)])
ns=M.nullspace(); print('R=',R,'kernel dim',len(ns))
v=sp.simplify(ns[0]); v=(v/ v[R.index((2,1))]).applyfunc(sp.factor); print('pi_q kernel:',list(v))
v1=v.applyfunc(lambda e: sp.limit(e,q,1)); print('q->1:',list(v1))
M1=sp.Matrix([[sum(chw(r,j,m).values()) for r in R] for j in range(10)]); n1=M1.nullspace()[0]; n1=n1/n1[R.index((2,1))]
print('x=1 kernel:',list(n1), 'equal:', list(v1)==list(n1))
# Hilbert polys at x=1
J=sp.Symbol('j'); print({r:sp.factor(sp.interpolate([sum(chw(r,j,m).values()) for j in range(8)],J)) for r in R})

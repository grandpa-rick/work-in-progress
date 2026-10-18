import sys, pickle, itertools; sys.path.insert(0,'.')
import sympy as sp
from fast import subsets, geq
N=int(sys.argv[1]); d=pickle.load(open('dual_N%d.pkl'%N,'rb'))
w=list(sp.symbols('w1:%d'%(N+1))); S=sorted(subsets(N),key=lambda T:(len(T),sum(T)))
comp=lambda M: tuple(c for c in range(1,N+1) if c not in M)
inv={wi:1/wi for wi in w}
for wxi in (0,1):
  for wy in (0,1):
    for k in range(-4,5):
      Y=sp.Matrix(len(S),len(S),lambda i,j: sp.cancel(d[wy,'A',comp(S[j]),comp(S[i])].subs(inv,simultaneous=True)*sp.prod([w[m-1]**k for m in S[j]])))
      X=sp.Matrix(len(S),len(S),lambda i,j: d[wxi,'xi',S[i],S[j]])
      C=(X*Y.inv()).applyfunc(sp.cancel)
      # allow column-scaling: normalise each row of C by diagonal
      ok=all(sp.cancel(C[i,j]/C[i,i]).is_number for i in range(len(S)) for j in range(len(S)) if C[i,j]!=0)
      if ok:
        print('xi',wxi,'Y',wy,'k',k,'C normalised constant!')
        for i in range(len(S)): print('  ',S[i],[sp.cancel(C[i,j]/C[i,i]) for j in range(len(S))], 'diag',sp.factor(C[i,i]))

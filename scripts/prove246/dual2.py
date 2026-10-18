import sys, itertools, pickle; sys.path.insert(0,'.')
import sympy as sp
from fast import subsets, geq
from rho_sym import Fs
N=int(sys.argv[1])
w=list(sp.symbols('w1:%d'%(N+1)))
S=sorted(subsets(N),key=lambda T:(len(T),sum(T)))
def eM(M):
    eps=[1 if i+1 in M else -1 for i in range(N)]
    return sp.prod([1-w[i]**(-eps[i])*w[j]**(-eps[j]) for i in range(N) for j in range(i,N)])
out={}
for wt in (0,1):
    A=sp.Matrix(len(S),len(S),lambda i,j: Fs(S[j],S[i],w,wt,N))
    Ainv=A.inv()
    for i,nu in enumerate(S):
        for j,M in enumerate(S):
            out[wt,'A',M,nu]=A[j,i]
            out[wt,'xi',nu,M]=sp.factor(sp.cancel(Ainv[i,j]*eM(M)))
pickle.dump(out,open('dual_N%d.pkl'%N,'wb'))
comp=lambda M: tuple(c for c in range(1,N+1) if c not in M)
inv={wi:1/wi for wi in w}
for wt in (0,1):
  for wt2 in (0,1):
    print('xi wt',wt,'vs A wt2',wt2,'at complement, w->1/w')
    for nu in S:
        rats=[]
        for M in S:
            X=out[wt,'xi',nu,M]; Y=out[wt2,'A',comp(M),comp(nu)].subs(inv,simultaneous=True)
            if X==0 and Y==0: continue
            rats.append(sp.factor(sp.cancel(X/Y)) if Y!=0 else 'Y0')
        print(' ',nu,rats)

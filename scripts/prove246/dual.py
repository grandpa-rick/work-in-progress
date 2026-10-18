import sys, itertools; sys.path.insert(0,'.')
import sympy as sp
from fast import subsets, geq
from rho_sym import Fs
N=int(sys.argv[1]); wt=int(sys.argv[2])
w=list(sp.symbols('w1:%d'%(N+1)))
S=sorted(subsets(N),key=lambda T:(len(T),sum(T)))
A=sp.Matrix(len(S),len(S),lambda i,j: Fs(S[j],S[i],w,wt,N))
Ainv=A.inv()  # triangular
def eM(M):
    eps=[1 if i+1 in M else -1 for i in range(N)]
    return sp.prod([1-w[i]**(-eps[i])*w[j]**(-eps[j]) for i in range(N) for j in range(i,N)])
for nu in S:
    i=S.index(nu)
    print(nu, [sp.factor(sp.cancel(Ainv[i,j]*eM(S[j]))) for j in range(len(S))])

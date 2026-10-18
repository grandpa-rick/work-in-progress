"""Check: l^F_nu(e_L) * F_L(L) * phi0_{L^c}(L^c) / phi0_{nu^c}(L^c) = kappa(nu,L) monomial; find it."""
import sys, itertools, random; sys.path.insert(0,'.')
import sympy as sp
from fast import subsets, geq
from rho_sym import Fs
N=int(sys.argv[1]); w=list(sp.symbols('w1:%d'%(N+1)))
S=sorted(subsets(N),key=lambda T:(len(T),sum(T)))
comp=lambda M: tuple(c for c in range(1,N+1) if c not in M)
F={(M,I):Fs(I,M,w,1,0) for M in S for I in S}   # wt1, NN=0 -> exponent 2r
P0={(M,I):Fs(I,M,w,0,0) for M in S for I in S}
A=sp.Matrix(len(S),len(S),lambda i,j: F[S[i],S[j]])
Ai=A.inv()
for i,nu in enumerate(S):
    row=[]
    for j,L in enumerate(S):
        if Ai[i,j]==0: continue
        row.append((L,sp.factor(sp.cancel(Ai[i,j]*F[L,L]*P0[comp(L),comp(L)]/P0[comp(L),comp(nu)]))))
    print(nu,row)

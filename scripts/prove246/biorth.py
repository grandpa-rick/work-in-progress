import sys; sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake246')
from common import *
N=int(sys.argv[1]); A=load(N); S=order(subsets(N)); dim=N*(N+1)//2
w=sp.symbols('w1:%d'%(N+1))
import random
R=random.Random(1); pt={wi:sp.Rational(R.randint(2,50),R.randint(2,50)) for wi in w}
comp=lambda M: tuple(c for c in range(1,N+1) if c not in M)
P0={k:(sp.cancel(A[k].subs(pt).subs(q,0)) if A[k]!=0 else 0) for k in A}
Pi={k:(lead_inf(A[k].subs(pt),dim-sum(k[1])) if A[k]!=0 else 0) for k in A}
Ainf=sp.Matrix(len(S),len(S),lambda i,j: Pi[S[i],S[j]])   # rows M, cols I
Linv=Ainf.inv()   # Linv[nu,L] = l^inf_nu(e_L)
bad=0
for i,nu in enumerate(S):
    for j,L in enumerate(S):
        pred=P0[comp(L),comp(nu)]/(Pi[L,L]*P0[comp(L),comp(L)])
        if sp.cancel(pred-Linv[i,j])!=0: bad+=1
print('N',N,'mismatches',bad,'of',len(S)**2)
# and the q=0 version
A0=sp.Matrix(len(S),len(S),lambda i,j: P0[S[i],S[j]]); L0=A0.inv(); bad=0
for i,nu in enumerate(S):
    for j,L in enumerate(S):
        pred=Pi[comp(L),comp(nu)]/(P0[L,L]*Pi[comp(L),comp(L)])
        if sp.cancel(pred-L0[i,j])!=0: bad+=1
print('q=0 version mismatches',bad)

# Independent check of Prop 2 at N=4: phi^inf from the stored full-q A (lead_inf), phi0 = A|_{q=0}; no closed form used.
import sys; sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake246')
from common import *
import random
N=4; A=load(N); S=order(subsets(N)); dim=N*(N+1)//2
w=sp.symbols('w1:%d'%(N+1)); R=random.Random(2)
pt={wi:sp.Rational(R.randint(2,40),R.randint(2,40)) for wi in w}
comp=lambda M: tuple(c for c in range(1,N+1) if c not in M)
Aq={k:(A[k].subs(pt) if A[k]!=0 else 0) for k in A}
P0={k:(sp.cancel(Aq[k].subs(q,0)) if Aq[k]!=0 else 0) for k in Aq}
Pi={k:(lead_inf(Aq[k],dim-sum(k[1])) if Aq[k]!=0 else 0) for k in Aq}
M=sp.Matrix(len(S),len(S),lambda i,j: Pi[S[i],S[j]]); Li=M.inv()
bad=sum(1 for i,nu in enumerate(S) for j,L in enumerate(S) if sp.cancel(P0[comp(L),comp(nu)]/(Pi[L,L]*P0[comp(L),comp(L)])-Li[i,j])!=0)
print('N=4 independent Prop 2 mismatches:',bad,'of',len(S)**2)

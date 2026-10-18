import sys, ast; sys.path.insert(0,'.'); sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake246')
from common import *; from cf import Phi
N=int(sys.argv[1]); A=load(N); S=subsets(N); w=sp.symbols('w1:%d'%(N+1)); dim=N*(N+1)//2
rat={}
for M in S:
  for I in S:
    if A[M,I]==0: assert Phi(I,M,w)==0; continue
    rat[M,I]=sp.cancel(lead_inf(A[M,I],dim-sum(I))/Phi(I,M,w))
# row*col test: rat[M,I]/rat[M,()] independent of M ; rat[M,()] = c*prod_{m in M} w_m^e
for I in S:
  print(I, set(sp.cancel(rat[M,I]/rat[M,()]) for M in S if (M,I) in rat))
print('row', [sp.factor(rat[M,()]/rat[(),()]) for M in S])
for e in range(-3,4):
  ok=True
  for I in S:
    vs=set(sp.cancel(rat[M,I]/rat[M,()]/sp.prod([w[m-1]**(e*len(I)) for m in M])) for M in S if (M,I) in rat)
    if len(vs)>1: ok=False
  print(e, ok)

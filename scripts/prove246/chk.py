import sys; sys.path.insert(0,'.'); sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake246')
from common import *; from cf import phiinf_cf
N=int(sys.argv[1]); A=load(N); S=subsets(N); w=sp.symbols('w1:%d'%(N+1)); dim=N*(N+1)//2
for M in S:
  for I in S:
    if A[M,I]==0: continue
    L=lead_inf(A[M,I],dim-sum(I)); C=phiinf_cf(I,M,w)
    print(M,I,sp.factor(sp.cancel(L/C)))

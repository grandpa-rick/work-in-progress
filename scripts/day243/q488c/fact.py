import sympy as sp,pickle,sys,itertools
from engine import subsets, geq
N=int(sys.argv[1]); A0,A1=pickle.load(open('A01_N%d.pkl'%N,'rb')); w=sp.symbols('w1:%d'%(N+1))
for a,b in itertools.combinations(range(1,N+1),2):
  low=[x for x in range(1,a)]
  S=[tuple(c) for k in range(len(low)+1) for c in itertools.combinations(low,k)]
  for Mpp in S:
    M=tuple(sorted(Mpp+(a,b)))
    rat=set()
    for kpp in S:
      if not geq(Mpp,kpp): continue
      K=tuple(sorted(kpp+(a,b)))
      rat.add(sp.factor(sp.cancel(A0[M,K]/A0[Mpp,kpp])))
    print((a,b),Mpp,len(rat),rat if len(rat)==1 else '')

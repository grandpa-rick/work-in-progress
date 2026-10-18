import sys; sys.path.insert(0,'.'); sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake246')
from common import *
from rho_sym import Fs
N=int(sys.argv[1]); A=load(N); S=subsets(N); dim=N*(N+1)//2
w=sp.symbols('w1:%d'%(N+1))
import random
R=random.Random(5); pt={wi:sp.Rational(R.randint(2,50),R.randint(2,50)) for wi in w}
def c(I):
    r=1
    for i in I: r*=(-1)**i*sp.prod(w[:i])
    return r
vals=set()
for M in S:
  for I in S:
    if A[M,I]==0: continue
    L=lead_inf(A[M,I].subs(pt),dim-sum(I))
    Fv=Fs(I,M,list(w),1,0).subs(pt)
    pred=sp.sympify(sp.prod([w[m-1]**(-2*N) for m in M])).subs(pt)*Fv/sp.sympify(c(I)).subs(pt)
    vals.add(sp.nsimplify(L/pred))
print(N, vals)

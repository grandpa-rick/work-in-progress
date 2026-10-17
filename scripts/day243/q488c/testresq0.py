import sympy as sp,pickle,itertools
from resq0 import residue
N=3; w=sp.symbols('w1:4'); res=pickle.load(open('res_N3.pkl','rb')); ok=bad=0
for (l,m,nu),(j,c0,c1,lp0,lp1) in res.items():
  for a,b in itertools.combinations(nu,2):
    h=1-1/(w[a-1]*w[b-1])
    R=sp.cancel(sp.cancel(c1*h).subs(w[b-1],1/w[a-1]))
    Rn=residue(l,m,nu,a,b,N,symbolic=True)
    x=sp.Symbol('x')
    if sp.cancel(R.subs(w[a-1],x)-Rn)==0: ok+=1
    else: bad+=1; print('MISMATCH',l,m,nu,(a,b),R,Rn)
print('ok',ok,'bad',bad)

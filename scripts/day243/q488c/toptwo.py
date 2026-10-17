import sympy as sp,pickle,sys,itertools
from engine import subsets, geq
from q0formula import phi0
N=int(sys.argv[1]); w=sp.symbols('w1:%d'%(N+1))
res=pickle.load(open('res_N%d.pkl'%N,'rb'))
def coord(F, nu, n):  # coordinates of function F on subsets of [n] in basis phi0, return coord at nu
    S=[s for s in subsets(n)] if n>0 else [()]
    order=sorted(S,key=lambda T:(len(T),sum(T))); c={}
    for x in order:
        c[x]=sp.cancel((F(x)-sum(c[k]*phi0(k,x,w) for k in c if c[k]!=0))/phi0(x,x,w))
        if x==nu: return c[x]
ok=bad=0
for (l,m,nu),(j,c0,c1,lp0,lp1) in res.items():
    if len(nu)<2: continue
    a,b=nu[-2],nu[-1]; nupp=nu[:-2]
    h=1-1/(w[a-1]*w[b-1])
    R=sp.cancel(sp.cancel(c1*h).subs(w[b-1],1/w[a-1]))
    if not (all(x<a for x in l) and all(x<a for x in m)): pred=0
    else:
        L=lambda M: sp.prod([w[c-1]**2 for c in M])
        F=lambda M: L(M)*phi0(l,M,w)*phi0(m,M,w)
        Xi=coord(F,nupp,a-1)
        K=sp.prod([(1-w[c-1]/w[a-1])*(1-w[a-1]*w[c-1]) for c in range(1,a)])*sp.prod([1-w[a-1]*w[c-1] for c in range(a+1,b)])
        pred=-Xi/K
    if sp.cancel(R-pred)==0: ok+=1
    else: bad+=1; print('MISMATCH',l,m,nu,sp.factor(R),sp.factor(pred))
print('top-two residue formula: ok',ok,'bad',bad)

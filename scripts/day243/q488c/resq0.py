"""Res_{H_ab} c1^nu_{lam,mu} via proved formula  Res = gamma * Res_H l^0_nu(1_{⊇ab} phi0_lam phi0_mu).
Uses only q=0 data on the lower ideal of nu. Other w's -> random rationals (seeded), w_a=x symbolic, w_b=1/x."""
import sympy as sp, itertools, random, sys
from engine import geq
from q0formula import phi0
def lower(nu,N):
    return [tuple(c) for k in range(len(nu)+1) for c in itertools.combinations(range(1,N+1),k) if geq(nu,tuple(c))]
def residue(lam,mu,nu,a,b,N,seed=1,symbolic=False):
    random.seed(seed)
    x,eps=sp.symbols('x eps')
    w=list(sp.symbols('w1:%d'%(N+1)))
    vals={}
    for c in range(1,N+1):
        if c not in (a,b) and not symbolic: vals[w[c-1]]=sp.Rational(random.randint(2,40),random.randint(41,90))*random.choice([1,-1])+random.randint(1,5)
    # approach H: w_a=x, w_b=(1+eps)/x? use h=1-1/(w_a w_b); set w_b = 1/(x*(1-eps)) so h=eps
    sub=dict(vals); sub[w[a-1]]=x; sub[w[b-1]]=1/(x*(1-eps))
    L=sorted(lower(nu,N),key=lambda T:(len(T),sum(T)))
    ph={}
    for k in L:
        for M in L:
            ph[k,M]=sp.cancel(phi0(k,M,w).subs(sub)) if geq(M,k) else 0
    F={M:(ph[lam,M] if lam in L else 0)*(ph[mu,M] if mu in L else 0)*(1 if (a in M and b in M) else 0) for M in L}
    c={}
    for M in L:
        c[M]=sp.cancel((F[M]-sum(c[k]*ph[k,M] for k in c if ph[k,M]!=0 and c[k]!=0))/ph[M,M])
    lnu=c[nu]
    r=sp.cancel(sp.cancel(lnu*eps).subs(eps,0))
    g=(x-1/x)**2
    return sp.factor(sp.cancel(g*r))
if __name__=='__main__':
    pass

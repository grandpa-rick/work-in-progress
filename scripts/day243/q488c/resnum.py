"""Numeric residue test: all w random rationals except eps (approach to H_ab: w_b = 1/(w_a(1-eps)), h=eps).
Res_H c1^nu = gamma * lim eps * l^0_nu(1_{⊇ab} phi0_lam phi0_mu)  (proved residue formula)."""
import sympy as sp, itertools, random
from fractions import Fraction as Fr
from engine import geq
def lower(nu,N):
    return [tuple(c) for k in range(len(nu)+1) for c in itertools.combinations(range(1,N+1),k) if geq(nu,tuple(c))]
eps=sp.Symbol('eps')
def phi0n(kappa,M,wv):
    r=len(kappa)
    if len(M)<r: return 0
    tot=0
    for h in itertools.permutations(M,r):
        hv=[wv[m] for m in h]; t=1
        for j in range(r):
            t*=(1-hv[j]**-2)
            for k in range(1,kappa[j]): t*=(1-wv[k]/hv[j])
            for k in range(j+1,r): t*=(1-1/(hv[j]*hv[k]))/(1-hv[j]/hv[k])
        for tm in M:
            if tm in h: continue
            for j in range(r): t*=(1-1/(wv[tm]*hv[j]))/(1-wv[tm]/hv[j])
        tot+=t
    return sp.cancel(tot)
def residue(lam,mu,nu,a,b,N,seed=1):
    random.seed(seed)
    wv={c:sp.Rational(random.randint(2,97),random.randint(2,97))+random.randint(0,3) for c in range(1,N+1)}
    wv[b]=1/(wv[a]*(1-eps))
    L=sorted(lower(nu,N),key=lambda T:(len(T),sum(T)))
    ph={(k,M):(phi0n(k,M,wv) if geq(M,k) else 0) for k in L for M in L}
    F={M:(ph[lam,M] if lam in L else 0)*(ph[mu,M] if mu in L else 0)*(1 if (a in M and b in M) else 0) for M in L}
    c={}
    for M in L:
        s=F[M]-sum(c[k]*ph[k,M] for k in c if ph[k,M]!=0 and c[k]!=0)
        c[M]=sp.cancel(s/ph[M,M])
    r=sp.cancel(c[nu]*eps).subs(eps,0)
    x=wv[a]; return sp.nsimplify((x-1/x)**2*r)

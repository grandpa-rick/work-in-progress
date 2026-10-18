"""closed-form candidates for phi^inf (Day 246 PROVE)."""
import sys, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake246')
def phiinf_cf(kappa, M, w):
    r=len(kappa); s=len(M); tot=0
    if s<r: return sp.Integer(0)
    for h in itertools.permutations(M, r):
        hv=[w[m-1] for m in h]; t=sp.Integer(1)
        for j in range(r):
            t*= hv[j]**(kappa[j]+s-r-1)*(1-hv[j]**-2)*sp.prod([1-w[k-1]/hv[j] for k in range(1,kappa[j])])
            for k in range(j+1,r): t*= (1-1/(hv[j]*hv[k]))/(1-hv[k]/hv[j])
        for tm in M:
            if tm in h: continue
            tv=w[tm-1]
            for j in range(r): t*=(1-1/(tv*hv[j]))/(1-hv[j]/tv)
        tot+=t
    tot*=(-1)**(r*s)*sp.prod([w[m-1]**r for m in M])
    return sp.factor(sp.cancel(tot))

def T0(kappa, h, M, w):
    """q=0 summand (2.1) for ordered head selection h (h[j] carries label kappa[j])."""
    r=len(kappa); hv=[w[m-1] for m in h]; t=sp.Integer(1)
    for j in range(r):
        t*= (1-hv[j]**-2)*sp.prod([1-w[k-1]/hv[j] for k in range(1,kappa[j])])
        for k in range(j+1,r): t*= (1-1/(hv[j]*hv[k]))/(1-hv[j]/hv[k])
    for tm in M:
        if tm in h: continue
        tv=w[tm-1]
        for j in range(r): t*=(1-1/(tv*hv[j]))/(1-tv/hv[j])
    return t
def Phi(kappa, M, w, wt=True):
    """wt=True: q=inf form sum_h T0 * prod h_j^{kappa_j-2j} (j 1-based); wt=False: phi0."""
    r=len(kappa)
    if len(M)<r: return sp.Integer(0)
    tot=0
    for h in itertools.permutations(M, r):
        t=T0(kappa,h,M,w)
        if wt: t*=sp.prod([w[h[j]-1]**(kappa[j]-2*(j+1)) for j in range(r)])
        tot+=t
    return sp.factor(sp.cancel(tot))

def F(kappa, M, w, wt=True):
    """gauge-normalised basis: wt=True -> phi^inf (up to row c*prod_{m in M}w_m^e and column monomial), wt=False -> phi0."""
    v=Phi(kappa,M,w,wt)
    if wt: v*=sp.prod([w[m-1]**(2*len(kappa)) for m in M])
    return v

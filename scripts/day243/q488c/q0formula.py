import sympy as sp, itertools
def phi0(kappa, M, w):
    r=len(kappa); tot=0
    if len(M)<r: return sp.Integer(0)
    for h in itertools.permutations(M, r):
        hv=[w[m-1] for m in h]; t=sp.Integer(1)
        for j in range(r):
            t*= (1-hv[j]**-2)*sp.prod([1-w[k-1]/hv[j] for k in range(1,kappa[j])])
            for k in range(j+1,r): t*= (1-1/(hv[j]*hv[k]))/(1-hv[j]/hv[k])
        for tm in M:
            if tm in h: continue
            tv=w[tm-1]
            for j in range(r): t*=(1-1/(tv*hv[j]))/(1-tv/hv[j])
        tot+=t
    return sp.factor(sp.cancel(tot))

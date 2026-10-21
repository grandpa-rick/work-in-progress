# Demazure on OG(n,2n+1), type B_n, doubled coordinates (weights in 1/2 Z^n -> Z^n)
import itertools, sympy as sp
from collections import defaultdict
from dem import word, chi as chiC
def DiB(f,i,n):
    g=defaultdict(int)
    for mu,a in f.items():
        mu=list(mu)
        if i<n: k=(mu[i-1]-mu[i])//2; al=[0]*n; al[i-1]=2; al[i]=-2
        else: k=mu[n-1]; al=[0]*n; al[n-1]=2   # <mu,2e_n> with mu doubled: mu_n(doubled)=2*mu_n -> 2mu_n = mu[n-1]
        if k>=0:
            for t in range(k+1): g[tuple(mu[j]-t*al[j] for j in range(n))]+=a
        elif k<=-2:
            for t in range(1,-k): g[tuple(mu[j]+t*al[j] for j in range(n))]-=a
    return {k:v for k,v in g.items() if v}
def chiB(kap,n,m):
    ws=word(kap,n); f={tuple([m]*n):1}   # m*omega_n^B = (m/2,..), doubled -> (m,..,m)
    for i in reversed(ws): f=DiB(f,i,n)
    return sum(f.values())
if __name__=="__main__":
    import sys
    n=int(sys.argv[1]); x=sp.Symbol('m'); c=sp.Symbol('c')
    for r in range(n+1):
        for kap in itertools.combinations(range(1,n+1),r):
            d=sum(kap)
            PC=sp.interpolate([(m,chiC(set(kap),n,m)[0]) for m in range(d+2)],x)
            PB=sp.interpolate([(m,chiB(set(kap),n,m)) for m in range(d+2)],x)
            pc=sp.factor((-1)**d*PC.subs(x,-c/2-1)); pb=sp.factor((-1)**d*PB.subs(x,-c-1))
            print(list(kap),'C:',pc,'  B:',pb,'  ratio',sp.simplify(pc/pb) if pb!=0 else None)

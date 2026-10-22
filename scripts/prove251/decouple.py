import sys, itertools, sympy as sp
from pq import polys, x as X
z=sp.Symbol('z')
def A(N,s,ks,k):  # [z^k] (1+z)^{N-s}(1-z)^{ks-N+s}
    if k<0: return 0
    return sp.series((1+z)**(N-s)*(1-z)**(ks-N+s),z,0,k+1).removeO().coeff(z,k)
def J(S,lam,N):
    # nested coefficient extraction [prod z_s^{kappa_s}] prod_s f_s(z_s) prod_{r<s in S} (z_r-z_s)/(z_r+z_s), |z_s|<|z_r| for r<s
    S=list(S)
    if not S: return 1
    # expand kernel: prod over pairs of (1+2 sum_{n>=1} (-1)^n (z_s/z_r)^n): enumerate shift vectors
    pairs=[(r,s) for i,r in enumerate(S) for s in S[i+1:]]
    tot=0
    maxn=sum(lam[s-1] for s in S)
    def rec(pi, shift, coef):
        nonlocal tot
        if pi==len(pairs):
            v=coef
            for s in S:
                v*=A(N,s,lam[s-1],lam[s-1]-shift[s]) * sp.Rational(1,2**(lam[s-1]+1))
                if v==0: return
            tot+=v; return
        r,s=pairs[pi]
        for n in range(0,maxn+1):
            sh=dict(shift); sh[r]=sh[r]-n; sh[s]=sh[s]+n   # z_s^n z_r^{-n}: need [z_r^{k_r}] of f_r * z_r^{-n} -> A_r(k_r+n); z_s: A_s(k_s-n)
            c=coef*(1 if n==0 else 2*(-1)**n)
            rec(pi+1,sh,c)
    rec(0,{s:0 for s in S},1)
    return tot
def Hdec(kap,c):
    lam=sorted(kap,reverse=True); l=len(lam); N=int(2*c)
    tot=0
    for S in itertools.chain.from_iterable(itertools.combinations(range(1,l+1),q) for q in range(l+1)):
        pre=sp.Mul(*[-sp.Rational((-1)**lam[r-1],2**(lam[r-1]+1)) for r in range(1,l+1) if r not in S])
        tot+=(-1)**len(S)*pre*J(S,lam,N)
    return (-1)**(sum(lam)-l)*tot
ok=tot=0
for kap in [(1,),(3,),(1,2),(1,3),(2,4),(1,2,3),(1,2,4),(2,3,4),(1,2,3,4)]:
    S_,E,O,p,q=polys(kap)
    for c in [sp.Rational(n,2) for n in range(1,9)]:
        v=Hdec(kap,c); w=E.subs(X,-2*c); tot+=1; ok+=(v==w)
        if v!=w: print('MISS',kap,c,v,w)
print('decoupled formula',ok,'/',tot)

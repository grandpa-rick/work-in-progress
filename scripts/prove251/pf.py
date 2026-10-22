import sys, itertools
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb
def gen_binom(e,k):  # [z^k](1+z)^e or (1-z)^e general integer e: coefficient of z^k in (1+z)^e
    if k<0: return 0
    num=1
    for i in range(k): num*=(e-i)
    from math import factorial
    return Fr(num,factorial(k))
def coeff_f(N,s,ks,k):  # [z^k](1+z)^{N-s}(1-z)^{ks-N+s}
    if k<0: return 0
    return sum(gen_binom(N-s,i)*gen_binom(ks-N+s,k-i)*(-1)**(k-i) for i in range(k+1))
def g_coeff(N,s,ks,k):  # [z^k] g_s, g_s = -[(-1)^ks z^ks + f_s]
    v=-coeff_f(N,s,ks,k)
    if k==ks: v-=(-1)**ks
    return v
def Hhat_direct(lam,N):
    l=len(lam); pairs=[(r,s) for r in range(l) for s in range(r+1,l)]
    maxn=sum(lam)
    tot=Fr(0)
    def rec(pi,shift,coef):
        nonlocal tot
        if pi==len(pairs):
            v=coef
            for s in range(l):
                v*=g_coeff(N,s+1,lam[s],lam[s]-shift[s])
                if v==0: return
            tot+=v; return
        r,s=pairs[pi]
        for n in range(0,maxn+1):
            sh=list(shift); sh[r]-=n; sh[s]+=n
            rec(pi+1,sh,coef*(1 if n==0 else 2*(-1)**n))
    rec(0,[0]*l,Fr(1)); return tot
def pair_entry(N,r,s,kr,ks):  # [z_r^kr z_s^ks] g_r g_s (z_r-z_s)/(z_r+z_s), |z_s|<|z_r|
    return sum((1 if n==0 else 2*(-1)**n)*g_coeff(N,r,kr,kr+n)*g_coeff(N,s,ks,ks-n) for n in range(0,ks+1))
def pfaff(M):
    n=len(M)
    if n==0: return Fr(1)
    if n%2: return Fr(0)
    tot=Fr(0)
    for j in range(1,n):
        if M[0][j]==0: continue
        idx=[k for k in range(1,n) if k!=j]
        sub=[[M[a][b] for b in idx] for a in idx]
        tot+=(-1)**(j-1)*M[0][j]*pfaff(sub)
    return tot
def Hhat_pf(lam,N):
    l=len(lam); idx=list(range(1,l+1)); n=l if l%2==0 else l+1
    M=[[Fr(0)]*n for _ in range(n)]
    for a in range(l):
        for b in range(a+1,l):
            M[a][b]=pair_entry(N,a+1,b+1,lam[a],lam[b]); M[b][a]=-M[a][b]
        if l%2:
            M[a][l]=g_coeff(N,a+1,lam[a],lam[a]); M[l][a]=-M[a][l]
    return pfaff(M)
if __name__=='__main__':
    sys.path.insert(0,'.')
    from pq import polys, x as X
    import sympy as sp
    ok1=ok2=tot=0
    for kap in [(1,),(3,),(1,2),(1,3),(2,4),(1,2,3),(1,2,4),(2,3,4),(1,2,3,4),(1,3,4),(2,3,5),(1,2,3,5)]:
        lam=sorted(kap,reverse=True); l=len(lam)
        S_,E,O,p,q=polys(kap)
        for N in range(1,13):
            w=E.subs(X,-N)*(-1)**(sum(lam)-l)*sp.Integer(2)**sum(k+1 for k in lam)
            d=Hhat_direct(lam,N) if sum(lam)<=8 else None
            pf=Hhat_pf(lam,N)
            tot+=1; ok2+= (sp.Rational(pf.numerator,pf.denominator)==w)
            if d is not None: ok1+= (sp.Rational(d.numerator,d.denominator)==w)
            else: ok1+=1
            if sp.Rational(pf.numerator,pf.denominator)!=w: print('PF MISS',kap,N,pf,w)
    print('direct',ok1,'/',tot,' pfaffian',ok2,'/',tot)

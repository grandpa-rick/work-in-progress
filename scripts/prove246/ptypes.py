"""Direct (dual-free) evaluation of the three problem types, for checking the reduction steps.
A: Res_{H_ab} l^F_nu(1_{>=ab} w^c F_lam F_mu)      (system [n], w_b = 1/(x(1-eps)))
B: l^F_kap(G), G(M)=[a in M] w_S^c F_lam F_mu(S)/g_a(S), S=M\\a, g_a(S)=prod s^2(1-x/s)/(1-x s), x=w_a
C: l^F_kap(w^c F_lam F_mu)"""
import sys, itertools, random; sys.path.insert(0,'.')
import fast
from fast import *
from fractions import Fraction as Fr
def mk(vals,a=None,b=None):
    w=[S.const(v) for v in vals]
    if b is not None: w[b-1]=S(0,[Fr(1)/vals[a-1]]*fast.K)
    return w
class Sys:
    def __init__(s,n,w): s.n=n; s.w=w; s.cache={}
    def F(s,I,M):
        if (M,I) not in s.cache: s.cache[M,I]=Fv(I,M,s.w,1,NN=0) if geq(M,I) else S.const(0)
        return s.cache[M,I]
    def coord(s,f,nu):
        down=sorted([M for M in subsets(s.n) if geq(nu,M)],key=lambda T:(len(T),sum(T)))
        c={}
        for M in down:
            v=f.get(M,S.const(0))
            for I in c:
                if geq(M,I) and not c[I].iszero(): v=v-c[I]*s.F(I,M)
            c[M]=v/s.F(M,M)
        return c[nu]
    def wpow(s,M,c):
        r=S.const(1)
        for m in M: r=r*(s.w[m-1]**c)
        return r
def A(n,vals,a,b,lam,mu,nu,c):
    s=Sys(n,mk(vals,a,b)); f={}
    for M in subsets(n):
        if a in M and b in M and geq(nu,M): f[M]=s.wpow(M,c)*s.F(lam,M)*s.F(mu,M)
    return s.coord(f,nu).coef(-1)
def B(n,vals,a,lam,mu,kap,c):
    s=Sys(n,mk(vals)); x=s.w[a-1]; f={}; one=S.const(1)
    for M in subsets(n):
        if a in M and geq(kap,M):
            Sx=tuple(m for m in M if m!=a); g=S.const(1)
            for m in Sx: g=g*(s.w[m-1]**2)*(one-x/s.w[m-1])/(one-x*s.w[m-1])
            f[M]=s.wpow(Sx,c)*s.F(lam,Sx)*s.F(mu,Sx)/g
    return s.coord(f,kap).coef(0)
def C(n,vals,lam,mu,kap,c):
    s=Sys(n,mk(vals)); f={}
    for M in subsets(n):
        if geq(kap,M): f[M]=s.wpow(M,c)*s.F(lam,M)*s.F(mu,M)
    return s.coord(f,kap).coef(0)

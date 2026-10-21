"""Modular (P=2^61-1) engine for A/B/C weighted residues. Series in eps truncated at K terms (relative)."""
import itertools, random
from fractions import Fraction as Fr
P=None
K=6
def inv(a): return 1/Fr(a)
class S:
    __slots__=('v','c')
    def __init__(s,v,c): s.v=v; s.c=c
    @staticmethod
    def const(a): return S(0,[Fr(a)]+[0]*(K-1))
    def norm(s):
        c=s.c; i=0
        while i<K and c[i]==0: i+=1
        if i==K: return S(0,[0]*K)
        return S(s.v+i,c[i:]+[0]*i)
    def iszero(s): return all(x==0 for x in s.c)
    def __add__(s,o):
        if not isinstance(o,S): o=S.const(o)
        if s.iszero(): return o
        if o.iszero(): return s
        v=min(s.v,o.v); c=[0]*K; ds=s.v-v; do=o.v-v
        for i in range(K):
            x=0
            if 0<=i-ds: x+=s.c[i-ds]
            if 0<=i-do: x+=o.c[i-do]
            c[i]=x
        return S(v,c).norm()
    __radd__=__add__
    def __neg__(s): return S(s.v,[(-x) for x in s.c])
    def __sub__(s,o): return s+(-o if isinstance(o,S) else S.const(-o))
    def __mul__(s,o):
        if not isinstance(o,S): return S(s.v,[x*o for x in s.c])
        c=[0]*K
        for i in range(K):
            if s.c[i]==0: continue
            for j in range(K-i): c[i+j]=(c[i+j]+s.c[i]*o.c[j])
        return S(s.v+o.v,c).norm()
    __rmul__=__mul__
    def inv(s):
        s=s.norm(); assert not s.iszero()
        c=[0]*K; i0=inv(s.c[0]); c[0]=i0
        for n in range(1,K):
            c[n]=(-sum(s.c[k]*c[n-k] for k in range(1,n+1))*i0)
        return S(-s.v,c)
    def __truediv__(s,o):
        if not isinstance(o,S): return s*inv(o)
        return s*o.inv()
    def __pow__(s,n):
        if n<0: return s.inv()**(-n)
        r=S.const(1); b=s
        while n:
            if n&1: r=r*b
            b=b*b; n>>=1
        return r
    def coef(s,k):
        i=k-s.v; return s.c[i] if 0<=i<K else 0
def subsets(N):
    out=[]
    for s in range(N+1): out+=[tuple(c) for c in itertools.combinations(range(1,N+1),s)]
    return out
def geq(M,I):
    s,r=len(M),len(I); return s>=r and all(M[s-r+j]>=I[j] for j in range(r))
def T0(kappa,h,M,w):
    r=len(kappa); hv=[w[m-1] for m in h]; t=S.const(1); one=S.const(1)
    for j in range(r):
        t=t*(one-hv[j]**-2)
        for k in range(1,kappa[j]): t=t*(one-w[k-1]/hv[j])
        for k in range(j+1,r): t=t*(one-(hv[j]*hv[k]).inv())/(one-hv[j]/hv[k])
    for tm in M:
        if tm in h: continue
        tv=w[tm-1]
        for j in range(r): t=t*(one-(tv*hv[j]).inv())/(one-tv/hv[j])
    return t
def Fv(kappa,M,w):
    r=len(kappa); tot=S.const(0)
    if len(M)<r: return tot
    for h in itertools.permutations(M,r):
        if any(h[j]<kappa[j] for j in range(r)): continue
        t=T0(kappa,h,M,w)
        for j in range(r): t=t*(w[h[j]-1]**(kappa[j]-2*(j+1)))
        tot=tot+t
    for m in M: tot=tot*(w[m-1]**(2*r))
    return tot
class Sys:
    def __init__(s,n,seed,a=None,b=None):
        R=random.Random(seed); s.n=n
        vals=[R.randrange(2,P-1) for _ in range(n)]
        s.w=[S.const(v) for v in vals]
        if b is not None: s.w[b-1]=S(0,[inv(vals[a-1])]*K)
        s.cache={}; s.subs=subsets(n)
        s.down={}
    def F(s,I,M):
        k=(M,I)
        if k not in s.cache: s.cache[k]=Fv(I,M,s.w) if geq(M,I) else S.const(0)
        return s.cache[k]
    def coord(s,f,nu):
        if nu not in s.down: s.down[nu]=sorted([M for M in s.subs if geq(nu,M)],key=lambda T:(len(T),sum(T)))
        c={}
        for M in s.down[nu]:
            v=f.get(M,S.const(0))
            for I in c:
                if geq(M,I) and not c[I].iszero(): v=v-c[I]*s.F(I,M)
            c[M]=v/s.F(M,M)
        return c[nu]
    def wpow(s,M,c):
        r=S.const(1)
        for m in M: r=r*(s.w[m-1]**c)
        return r
def A(s,a,b,lam,mu,nu,c):
    f={}
    for M in s.subs:
        if a in M and b in M and geq(nu,M):
            Mpp=tuple(m for m in M if m not in (a,b))
            f[M]=s.wpow(Mpp,c)*s.F(lam,Mpp)*s.F(mu,Mpp)
    return s.coord(f,nu).coef(-1)
def B(s,a,lam,mu,kap,c):
    x=s.w[a-1]; f={}; one=S.const(1)
    for M in s.subs:
        if a in M and geq(kap,M):
            Sx=tuple(m for m in M if m!=a); g=S.const(1)
            for m in Sx: g=g*(s.w[m-1]**2)*(one-x/s.w[m-1])/(one-x*s.w[m-1])
            f[M]=s.wpow(Sx,c)*s.F(lam,Sx)*s.F(mu,Sx)/g
    return s.coord(f,kap).coef(0)
def C(s,lam,mu,kap,c):
    f={}
    for M in s.subs:
        if geq(kap,M): f[M]=s.wpow(M,c)*s.F(lam,M)*s.F(mu,M)
    return s.coord(f,kap).coef(0)
def EA(n,a,b,lam,mu,nu):
    rest=[x for x in range(1,n+1) if x not in (a,b)]
    return [S_ for k in range(len(rest)+1) for S_ in itertools.combinations(rest,k) if geq(nu,tuple(sorted(S_+(a,b)))) and geq(S_,lam) and geq(S_,mu)]
def EB(n,a,lam,mu,kap):
    rest=[x for x in range(1,n+1) if x!=a]
    return [S_ for k in range(len(rest)+1) for S_ in itertools.combinations(rest,k) if geq(kap,tuple(sorted(S_+(a,)))) and geq(S_,lam) and geq(S_,mu)]
def EC(n,lam,mu,kap):
    return [M for M in subsets(n) if geq(kap,M) and geq(M,lam) and geq(M,mu)]

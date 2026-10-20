"""Exact A-residue engine over Q(w) (flint). Logic copied verbatim from prove247/sy.py / mp.py; only the field changed.
Variables: w1..wN, with w_b := 1/(w_a (1-eps)); residue = eps^{-1} coefficient. Result is in Q(w_j : j != b)."""
import itertools
from ff import Ctx, Fr
K=6
def setup(N):
    C=Ctx(['w%d'%i for i in range(1,N+1)]); Fr.C=C; return C
ONE=lambda: Fr(Fr.C.one_p)
ZERO=lambda: Fr(Fr.C.zero_p)
class S:
    __slots__=('v','c')
    def __init__(s,v,c): s.v=v; s.c=c
    @staticmethod
    def const(a):
        if not isinstance(a,Fr): a=Fr(Fr.C.one_p*a)
        return S(0,[a]+[ZERO() for _ in range(K-1)])
    def norm(s):
        c=s.c; i=0
        while i<K and c[i].n==0: i+=1
        if i==K: return S(0,[ZERO() for _ in range(K)])
        return S(s.v+i,c[i:]+[ZERO() for _ in range(i)])
    def iszero(s): return all(x.n==0 for x in s.c)
    def __add__(s,o):
        if not isinstance(o,S): o=S.const(o)
        if s.iszero(): return o
        if o.iszero(): return s
        v=min(s.v,o.v); c=[]; ds=s.v-v; do=o.v-v
        for i in range(K):
            x=ZERO()
            if 0<=i-ds: x=x+s.c[i-ds]
            if 0<=i-do: x=x+o.c[i-do]
            c.append(x)
        return S(v,c).norm()
    __radd__=__add__
    def __neg__(s): return S(s.v,[-x for x in s.c])
    def __sub__(s,o): return s+(-o if isinstance(o,S) else S.const(-o))
    def __mul__(s,o):
        if not isinstance(o,S): return S(s.v,[x*o for x in s.c])
        c=[ZERO() for _ in range(K)]
        for i in range(K):
            if s.c[i].n==0: continue
            for j in range(K-i):
                if o.c[j].n!=0: c[i+j]=c[i+j]+s.c[i]*o.c[j]
        return S(s.v+o.v,c).norm()
    __rmul__=__mul__
    def inv(s):
        s=s.norm(); assert not s.iszero()
        c=[ZERO() for _ in range(K)]; i0=s.c[0].inv(); c[0]=i0
        for n in range(1,K):
            acc=ZERO()
            for k in range(1,n+1):
                if s.c[k].n!=0: acc=acc+s.c[k]*c[n-k]
            c[n]=-(acc*i0)
        return S(-s.v,c)
    def __truediv__(s,o):
        if not isinstance(o,S): return s*(ONE()/o)
        return s*o.inv()
    def __pow__(s,n):
        if n<0: return s.inv()**(-n)
        r=S.const(1); b=s
        while n:
            if n&1: r=r*b
            b=b*b; n>>=1
        return r
    def coef(s,k):
        i=k-s.v; return s.c[i] if 0<=i<K else ZERO()
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
    def __init__(s,n,a,b):
        s.n=n; g=Fr.C.gens
        vals=[Fr(g[i]) for i in range(n)]
        s.w=[S.const(v) for v in vals]
        s.w[b-1]=S(0,[vals[a-1].inv() for _ in range(K)])
        s.cache={}; s.subs=subsets(n); s.down={}
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
def EA(n,a,b,lam,mu,nu):
    rest=[x for x in range(1,n+1) if x not in (a,b)]
    return [S_ for k in range(len(rest)+1) for S_ in itertools.combinations(rest,k) if geq(nu,tuple(sorted(S_+(a,b)))) and geq(S_,lam) and geq(S_,mu)]

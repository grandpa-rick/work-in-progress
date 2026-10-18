"""Fast exact engine: Laurent series in eps (w_b = 1/(x(1-eps))) over Fractions.
F_kappa(M): wt=0 -> phi0 (2.1); wt=1 -> phi^inf closed form  prod_M m^{2r} * sum_h T0(h) prod_j h_j^{kappa_j-2j}."""
import itertools, random
from fractions import Fraction as Fr
K=12
class S:
    __slots__=('v','c')
    def __init__(s,v,c): s.v=v; s.c=c
    @staticmethod
    def const(a): return S(0,[Fr(a)]+[Fr(0)]*(K-1))
    def norm(s):
        c=s.c; i=0
        while i<len(c) and c[i]==0: i+=1
        if i==len(c): return S(0,[Fr(0)]*K)
        return S(s.v+i, c[i:]+[Fr(0)]*i)   # NOTE: loses precision at tail; K large enough
    def iszero(s): return all(x==0 for x in s.c)
    def __add__(s,o):
        if not isinstance(o,S): o=S.const(o)
        if s.iszero(): return o
        if o.iszero(): return s
        v=min(s.v,o.v); c=[Fr(0)]*K
        for i in range(K):
            if 0<=i-(s.v-v)<K: c[i]+=s.c[i-(s.v-v)]
            if 0<=i-(o.v-v)<K: c[i]+=o.c[i-(o.v-v)]
        return S(v,c).norm()
    __radd__=__add__
    def __neg__(s): return S(s.v,[-x for x in s.c])
    def __sub__(s,o): return s+(-o if isinstance(o,S) else S.const(-o))
    def __mul__(s,o):
        if not isinstance(o,S): return S(s.v,[x*o for x in s.c])
        c=[Fr(0)]*K
        for i in range(K):
            if s.c[i]==0: continue
            for j in range(K-i): c[i+j]+=s.c[i]*o.c[j]
        return S(s.v+o.v,c).norm()
    __rmul__=__mul__
    def inv(s):
        s=s.norm(); assert not s.iszero()
        c=[Fr(0)]*K; c[0]=1/s.c[0]
        for n in range(1,K):
            c[n]=-sum(s.c[k]*c[n-k] for k in range(1,n+1))/s.c[0]
        return S(-s.v,c)
    def __truediv__(s,o):
        if not isinstance(o,S): return s*(Fr(1)/Fr(o))
        return s*o.inv()
    def __pow__(s,n):
        if n<0: return s.inv()**(-n)
        r=S.const(1)
        for _ in range(n): r=r*s
        return r
    def coef(s,k):  # coefficient of eps^k
        i=k-s.v; return s.c[i] if 0<=i<K else (Fr(0) if i<0 else None)
def point(N,a,b,seed):
    R=random.Random(seed)
    vals=[Fr(R.randint(2,97),R.randint(2,97))+R.randint(0,3) for _ in range(N)]
    w=[S.const(v) for v in vals]
    x=vals[a-1]
    w[b-1]=S(0,[Fr(1)/x]*K)   # 1/(x(1-eps)) = (1/x) sum eps^k
    return w, x
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
def Fv(kappa,M,w,wt,NN=None):
    NN=len(w) if NN is None else NN
    r=len(kappa); tot=S.const(0)
    if len(M)<r: return tot
    for h in itertools.permutations(M,r):
        if any(h[j]<kappa[j] for j in range(r)): continue   # identically zero terms
        t=T0(kappa,h,M,w)
        if wt:
            for j in range(r): t=t*(w[h[j]-1]**(kappa[j]-2*(j+1)))
        tot=tot+t
    if wt:
        for m in M: tot=tot*(w[m-1]**(2*r-2*NN))
    return tot
def subsets(N):
    out=[]
    for s in range(N+1): out+= [tuple(c) for c in itertools.combinations(range(1,N+1),s)]
    return out
def geq(M,I):
    s,r=len(M),len(I); return s>=r and all(M[s-r+j]>=I[j] for j in range(r))
class Basis:
    def __init__(s,N,a,b,seed,wt):
        s.N=N; s.w,s.x=point(N,a,b,seed); s.wt=wt; s.cache={}; s.a=a; s.b=b
    def A(s,M,I):
        if (M,I) not in s.cache: s.cache[M,I]=Fv(I,M,s.w,s.wt) if geq(M,I) else S.const(0)
        return s.cache[M,I]
    def coords(s,f,nu):
        """nu-coordinate of f (dict M->S), using only down-set of nu."""
        down=sorted([M for M in subsets(s.N) if geq(nu,M)],key=lambda T:(len(T),sum(T)))
        c={}
        for M in down:
            v=f.get(M,S.const(0))
            for I in c:
                if geq(M,I) and not c[I].iszero(): v=v-c[I]*s.A(M,I)
            c[M]=v/s.A(M,M)
        return c[nu]
    def resid(s,lam,mu,nu):
        a,b=s.a,s.b
        f={}
        for M in subsets(s.N):
            if a in M and b in M and geq(nu,M):
                f[M]=s.A(M,lam)*s.A(M,mu)
        return s.coords(f,nu).coef(-1)

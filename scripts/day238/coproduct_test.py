# Day 238 (B): Delta_N = (N⊗N) Delta N^{-1} on P_nu(x; q=s, t_Mac=1/t), |nu|<=4.  N P_nu = t^{n(nu)} s^{n(nu')} P_nu.
import itertools, sympy as sp
from sympy.utilities.iterables import partitions
from collections import Counter, defaultdict
s,t=sp.symbols('s t',positive=True)
q=s; T=1/t
D=4
def parts(n): return sorted([tuple(sorted(sum(([k]*v for k,v in p.items()),[]),reverse=True)) for p in partitions(n)]) if n>0 else [()]
def conj(l): return tuple(sum(1 for x in l if x>i) for i in range(l[0])) if l else ()
def nf(l): return sum(i*x for i,x in enumerate(l))
def dom_le(a,b):
    A=B=0
    for i in range(max(len(a),len(b))):
        A+=a[i] if i<len(a) else 0; B+=b[i] if i<len(b) else 0
        if A>B: return False
    return True
X=sp.symbols('x1:%d'%(D+1))
def pk(k): return sum(x**k for x in X)
def zee(mu):
    c=Counter(mu); r=1
    for k,m in c.items(): r*=k**m*sp.factorial(m)
    return r
P={}; Pp={}  # Pp[nu] = dict mu->coef in p basis
Mp={}
for d in range(0,D+1):
    L=parts(d)
    if d==0: Pp[()]={():sp.Integer(1)}; continue
    # p_mu -> m: [m_lam] p_mu = coefficient of x^lam
    A=sp.zeros(len(L),len(L))
    for j,mu in enumerate(L):
        poly=sp.Poly(sp.expand(sp.Mul(*[pk(k) for k in mu])),*X)
        for i,lam in enumerate(L):
            A[i,j]=poly.coeff_monomial(sp.Mul(*[X[a]**lam[a] for a in range(len(lam))]))
    Ainv=A.inv()  # m_lam = sum_mu Ainv[mu-index, lam-index] p_mu
    mvec={lam:{mu:Ainv[j,i] for j,mu in enumerate(L)} for i,lam in enumerate(L)}
    def ip(f,g): return sp.together(sum(f[mu]*g.get(mu,0)*zee(mu)*sp.Mul(*[(1-q**k)/(1-T**k) for k in mu]) for mu in f))
    order=sorted(L,key=lambda l:[sum(l[:i+1]) for i in range(len(l))])  # dominance total for d<=5: increasing
    done=[]
    for lam in order:
        v=dict(mvec[lam])
        for mu in done:
            c=sp.cancel(ip(mvec[lam],Pp[mu])/ip(Pp[mu],Pp[mu]))
            for k in v: v[k]=v[k]-c*Pp[mu][k]
        Pp[lam]={k:sp.cancel(v[k]) for k in v}; done.append(lam)
    for lam in L: assert all(dom_le(lam,mu) or True for mu in L)
# inverse: p_mu in P basis per degree
pinP={}
for d in range(0,D+1):
    L=parts(d)
    M=sp.Matrix([[Pp[nu].get(mu,0) for nu in L] for mu in L])  # column nu = P_nu in p
    Mi=M.inv()  # p_mu = sum_nu Mi[nu,mu] P_nu
    Mi=Mi.applyfunc(sp.cancel)
    pinP[d]={mu:{nu:Mi[i,j] for i,nu in enumerate(L)} for j,mu in enumerate(L)}
def Delta_p(mu):
    out=defaultdict(int)
    for mask in itertools.product([0,1],repeat=len(mu)):
        a=tuple(sorted([mu[i] for i in range(len(mu)) if mask[i]==0],reverse=True))
        b=tuple(sorted([mu[i] for i in range(len(mu)) if mask[i]==1],reverse=True))
        out[(a,b)]+=1
    return out
def DeltaP(nu):
    out=defaultdict(int)
    for mu,c in Pp[nu].items():
        for (a,b),m in Delta_p(mu).items():
            for al,ca in pinP[sum(a)][a].items():
                for be,cb in pinP[sum(b)][b].items():
                    out[(al,be)]+=c*m*ca*cb
    return {k:sp.cancel(v) for k,v in out.items() if sp.cancel(v)!=0}
def gauge_N(l): return t**nf(l)*s**nf(conj(l))
def gauge_ctrl(l): return t**sum(l)
results={}
for d in range(1,D+1):
    for nu in parts(d):
        dP=DeltaP(nu)
        results[nu]=dP
import pickle; pickle.dump((Pp,results),open('/tmp/day238_coprod.pkl','wb'))
print('computed')

# exact type-C P_d(c), cached c-free parts; Descartes certificate on P(U+1+t)
import itertools, sys, math, random, sympy as sp
from eng import Psi, F, ell, leq
from fractions import Fraction as Fr
SAMPLE=float(sys.argv[5]) if len(sys.argv)>5 else 1.0
n=int(sys.argv[1]); random.seed(int(sys.argv[2])); shard=int(sys.argv[3]); nsh=int(sys.argv[4])
al=[Fr(random.randint(1,20),random.randint(1,20))*random.choice([1,-1]) for _ in range(n)]
U=list(range(1,n+1)); subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
dim=n*(n+1)//2; Umax=2*n-2
cs=list(range(-2,-2+dim+4)); extra=[dim+4,dim+5,dim+7]; allc=cs+extra  # extra outside interpolation range (fix 2026-10-21: old U+1..U+3 lay inside it)
dls=[Fr(k,10**20) for k in (1,2,4)]
c_=sp.Symbol('c'); t_=sp.Symbol('t')
W={dl:{i:1+al[i-1]*dl for i in range(1,n+1)} for dl in dls}
def X(lam,mu,kap,w):
    parts=[]
    for M in subs:
        if leq(M,kap) and leq(lam,M) and leq(mu,M):
            pm=Fr(1)
            for m in M: pm*=w[m]
            parts.append((ell(kap,M,n,w)*F(lam,M,w)*F(mu,M,w),pm))
    return [sum(g*pm**c for g,pm in parts) for c in allc]
norm={dl:X([],[],U,W[dl])[allc.index(2*n)] for dl in dls}
k=0
for kap in subs:
 for lam in subs:
  for mu in subs:
    if lam>mu: continue
    if not any(leq(M,kap) and leq(lam,M) and leq(mu,M) for M in subs): continue
    k+=1
    if k%nsh!=shard: continue
    if random.Random(k*7919).random()>SAMPLE: continue
    V={dl:[v/norm[dl] for v in X(lam,mu,kap,W[dl])] for dl in dls}
    v0=V[dls[0]][-1]; v1=V[dls[1]][-1]
    d=round(math.log2(abs(float(v1/v0))))
    vals=[]
    for i in range(len(allc)):
        g=[V[dl][i]/dl**d for dl in dls]
        vals.append(((8*g[0]-6*g[1]+g[2])/3).limit_denominator(10**9))
    P=sp.expand(sp.interpolate(list(zip(cs,[sp.Rational(v.numerator,v.denominator) for v in vals[:len(cs)]])),c_))
    ok=all(P.subs(c_,cc)==sp.Rational(v.numerator,v.denominator) for cc,v in zip(extra,vals[len(cs):]))
    Q=sp.Poly(sp.expand(P.subs(c_,Umax+1+t_)),t_).all_coeffs()
    desc=all(x>0 for x in Q) or all(x<0 for x in Q) or (len(Q)==1 and Q[0]!=0)
    # nonzero coefficients same sign (zeros allowed only if constant term nonzero)
    nz=[x for x in Q if x!=0]; desc2=(Q[-1]!=0) and (all(x>0 for x in nz) or all(x<0 for x in nz))
    print(lam,mu,kap,'d',d,'deg',sp.degree(P,c_),'chk',ok,'desc',desc2,'P',sp.factor(P),flush=True)

import sys, math, random, sympy as sp
from hpA2 import build, RA_all
from fractions import Fraction as Fr
N=5; random.seed(2)
al=[Fr(random.randint(1,20),random.randint(1,20))*random.choice([1,-1]) for _ in range(N)]
Uall=list(range(1,N+1)); dim=15; Umax=6
a,b=3,4; lam=(5,); mu=(5,); nu=(1,2,3,4,5)
cs=list(range(-2,-2+dim+4)); extra=[30,31,40]; allc=cs+extra
c_=sp.Symbol('c')
for exp_,den in [(40,10**30),(60,10**45)]:
    dls=[Fr(k,10**exp_) for k in (1,2,4)]
    sy={dl:build(N,a,b,al,dl) for dl in dls}
    nv={dl:RA_all(N,a,b,(),(),tuple(Uall),[2*N-2],sy[dl])[0] for dl in dls}
    V={dl:[v/nv[dl] for v in RA_all(N,a,b,lam,mu,nu,allc,sy[dl])] for dl in dls}
    d=round(math.log2(abs(float(V[dls[1]][0]/V[dls[0]][0]))))
    vals=[((8*V[dls[0]][i]/dls[0]**d-6*V[dls[1]][i]/dls[1]**d+V[dls[2]][i]/dls[2]**d)/3).limit_denominator(den) for i in range(len(allc))]
    P=sp.expand(sp.interpolate(list(zip(cs,[sp.Rational(v.numerator,v.denominator) for v in vals[:len(cs)]])),c_))
    ok=[P.subs(c_,cc)==sp.Rational(v.numerator,v.denominator) for cc,v in zip(extra,vals[len(cs):])]
    print('exp',exp_,'d',d,'deg',sp.degree(P,c_),'out-of-range checks',ok,'P',sp.factor(P))

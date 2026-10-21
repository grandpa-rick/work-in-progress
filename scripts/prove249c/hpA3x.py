# exact version: P_d(c) via rationalized Richardson values + sympy interpolation
import itertools, sys, math, random, sympy as sp
from hpA2 import build, RA_all
from mpq import subsets, EA
from fractions import Fraction as Fr
N=int(sys.argv[1]); random.seed(int(sys.argv[2])); shard=int(sys.argv[3]); nsh=int(sys.argv[4])
al=[Fr(random.randint(1,20),random.randint(1,20))*random.choice([1,-1]) for _ in range(N)]
Uall=list(range(1,N+1)); subs=subsets(N); dim=N*(N+1)//2; Umax=2*N-4
cs=list(range(-2,-2+dim+4)); extra=[dim+4,dim+5,dim+7]   # outside interpolation range (fix 2026-10-21: old U+1..U+3 lay inside it)
dls=[Fr(k,10**20) for k in (1,2,4)]
c_=sp.Symbol('c')
import os, re
DONE=set()
if len(sys.argv)>5 and os.path.exists(sys.argv[5]):
    for line in open(sys.argv[5]):
        m=re.match(r'(\(\d+, \d+\)) (\(.*?\)) (\(.*?\)) (\(.*?\)) d ',line)
        if m: DONE.add(tuple(eval(m.group(i)) for i in (1,2,3,4)))
print('resuming, done',len(DONE),file=sys.stderr)
pairs=[(a,b) for a in Uall for b in Uall if a<b]
k=0
for (a,b) in pairs:
  sy=None
  for nu in subs:
    for lam in subs:
      for mu in subs:
        if lam>mu: continue
        if not EA(N,a,b,lam,mu,nu): continue
        k+=1
        if k%nsh!=shard: continue
        if ((a,b),lam,mu,nu) in DONE: continue
        if sy is None:
            sy={dl:build(N,a,b,al,dl) for dl in dls}
            normv={dl:RA_all(N,a,b,(),(),tuple(Uall),[2*N-2],sy[dl])[0] for dl in dls}
        allc=cs+extra
        V={dl:[v/normv[dl] for v in RA_all(N,a,b,lam,mu,nu,allc,sy[dl])] for dl in dls}
        i0=allc.index(Umax+1)
        v0=V[dls[0]][-1]; v1=V[dls[1]][-1]
        if v0==0: print('ZERO at U+3',(a,b),lam,mu,nu,flush=True); continue
        d=round(math.log2(abs(float(v1/v0))))
        vals=[]
        for i in range(len(allc)):
            g=[V[dl][i]/dl**d for dl in dls]
            vals.append(((8*g[0]-6*g[1]+g[2])/3).limit_denominator(10**9))
        P=sp.expand(sp.interpolate(list(zip(cs,[sp.Rational(v.numerator,v.denominator) for v in vals[:len(cs)]])),c_))
        ok=all(P.subs(c_,cc)==sp.Rational(v.numerator,v.denominator) for cc,v in zip(extra,vals[len(cs):]))
        deg=sp.degree(P,c_) if P!=0 else -1
        if deg>=1:
            try:
                fl=sp.factor_list(P)[1]; rts=[x for f,_ in fl if sp.degree(f,c_)>0 for x in sp.Poly(f,c_).nroots(n=20,maxsteps=500)]
                mr=max(float(sp.re(r)) for r in rts)
            except Exception: mr=float('nan')
            intz=[x for x in range(Umax+1,Umax+60) if P.subs(c_,x)==0]
        else: mr=None; intz=[]
        print((a,b),lam,mu,nu,'d',d,'deg',deg,'chk',ok,'maxRe',None if mr is None else round(mr,4),'intzeros>U',intz,'P',sp.factor(P) if deg>=1 else P,flush=True)

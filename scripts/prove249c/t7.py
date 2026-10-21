import sympy as sp, itertools, sys, random
from collections import defaultdict
import engs
from engs import *
n=int(sys.argv[1]); cs=eval(sys.argv[2]); ntr=int(sys.argv[3])
W=sp.symbols('w1:%d'%(n+1)); w={k:W[k-1] for k in range(1,n+1)}
U=list(range(1,n+1))
subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
cache={}
def lt(expr,xi):
    num,den=sp.fraction(sp.cancel(sp.together(expr)))
    def L(p):
        P=sp.Poly(sp.expand(p),*W); ts=P.terms()
        return max(ts,key=lambda t:sum(xi[i]*t[0][i] for i in range(n)))
    (en,cn),(ed,cd)=L(num),L(den)
    return tuple(a-b for a,b in zip(en,ed)), sp.Rational(cn,1)/cd
def TL(lam,mu,kap,M):
    key=(tuple(lam),tuple(mu),tuple(kap),tuple(M))
    if key not in cache:
        cache[key]=sp.cancel(sp.together(ell(kap,M,n,w)*F(lam,M,w)*F(mu,M,w)))
    return cache[key]
random.seed(3)
stats=defaultdict(lambda: defaultdict(int)); ex={}
for kap in subs:
  for lam in subs:
    for mu in subs:
      if lam>mu: continue
      E=[M for M in subs if leq(M,kap) and leq(lam,M) and leq(mu,M)]
      if not E: continue
      for t in range(ntr):
        xi=[random.gauss(0,1) for _ in range(n)]
        leads=[(lt(TL(lam,mu,kap,M),xi),M) for M in E]
        for c in cs:
          g=defaultdict(int); who=defaultdict(list)
          for ((e,co),M) in leads:
            e2=tuple(e[i]+(c if (i+1) in M else 0) for i in range(n)); g[e2]+=co; who[e2].append(M)
          top=max(g,key=lambda e:sum(xi[i]*e[i] for i in range(n)))
          k=(len(who[top]),g[top]); stats[c][k]+=1
          if len(who[top])>1 or abs(g[top])!=1: ex.setdefault((c,k),(lam,mu,kap,who[top]))
for c in cs: print(c,dict(stats[c]))
for k,v in list(ex.items())[:20]: print(k,v)

from cfg import *
import sys
n=int(sys.argv[1]); cs=eval(sys.argv[2]); ntr=int(sys.argv[3])
random.seed(7)
U=list(range(1,n+1))
subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
stats=defaultdict(lambda: defaultdict(int))
for kap in subs:
  for lam in subs:
    for mu in subs:
      if lam>mu: continue
      if not any(leq(M,kap) and leq(lam,M) and leq(mu,M) for M in subs): continue
      for c in cs:
        for t in range(ntr):
          xi=[0]+[random.gauss(0,1) for _ in range(n)]
          T=top(Cconfigs(lam,mu,kap,c,n,xi),xi)
          e,sig,cnt,Ls=T[0]
          stats[c][(cnt,sig)]+=1
for c in cs: print(c,dict(stats[c]))

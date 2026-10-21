from deg import *
import sys
from collections import defaultdict
n=int(sys.argv[1]); c=int(sys.argv[2]); nrand=int(sys.argv[3])
random.seed(11)
U=list(range(1,n+1))
subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
a=[Fr(random.randint(2,9),random.randint(2,9)) for _ in range(n)]
B=4*n+4*c
lexes=[]
for s in itertools.product([1,-1],repeat=n):
    lexes.append(('lex',s,[s[i]*B**i for i in range(n)]))   # w_n dominant
dirs=lexes+[('rand',None,[random.randint(-7,7) for _ in range(n)]) for _ in range(nrand)]
tot=0; nolex=0; none_=0; bad=[]
for kap in subs:
  for lam in subs:
    for mu in subs:
      if lam>mu: continue
      E=[M for M in subs if leq(M,kap) and leq(lam,M) and leq(mu,M)]
      if not E: continue
      tot+=1; oklex=False; okany=False; winners=set()
      for kind,s,xi in dirs:
        if kind=='rand' and okany: break
        sc=[]
        for M in E:
            d=deg_eval(hL(lam,mu,kap,M,n),xi,a)
            sc.append((d+c*sum(xi[i-1] for i in M),tuple(M)))
        sc.sort(reverse=True)
        if len(sc)==1 or sc[0][0]>sc[1][0]:
            okany=True; winners.add((kind,s,sc[0][1]))
            if kind=='lex': oklex=True
      if not oklex: nolex+=1
      if not okany: none_+=1; bad.append((lam,mu,kap))
print('n',n,'c',c,'instances',tot,'no lex dir works',nolex,'no dir works',none_)
print(bad[:20])

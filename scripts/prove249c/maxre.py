import sys, sympy as sp
from collections import Counter
c=sp.Symbol('c'); U=int(sys.argv[1]); cnt=Counter(); edge=[]; worst=-99
for fn in sys.argv[2:]:
  for line in open(fn):
    if ' P ' not in line: continue
    P=sp.sympify(line.split(' P ',1)[1])
    if not P.has(c): cnt['const']+=1; continue
    fl=sp.factor_list(P)[1]; r=[x for f,_ in fl if sp.degree(f,c)>0 for x in sp.Poly(f,c).nroots(n=20,maxsteps=500)]; m=max(float(sp.re(x)) for x in r)
    worst=max(worst,m)
    if abs(m-U)<1e-9: cnt['=U']+=1; edge.append(line.split(' d ')[0])
    elif m<U: cnt['<U']+=1
    else: cnt['>U']+=1; print('VIOL',line.strip())
print(dict(cnt),'max Re root overall',worst); print('edge instances (root exactly U):',edge[:12])

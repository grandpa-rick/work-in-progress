import sys, sympy as sp
from collections import Counter
c=sp.Symbol('c'); z=sp.Symbol('z')
U=int(sys.argv[1]); offs=eval(sys.argv[2]); fns=sys.argv[3:]
res={o:Counter() for o in offs}; ex={}
for fn in fns:
  for line in open(fn):
    if ' P ' not in line: continue
    P=sp.sympify(line.split(' P ',1)[1]); D=sp.degree(P,c) if P.has(c) else 0
    if D<=0: continue
    for o in offs:
      c0=o
      # h(z) = (1-z)^{D+1} * sum_k P(c0+2k) z^k, truncated
      ser=sum(P.subs(c,c0+2*k)*z**k for k in range(D+3))
      h=sp.Poly(sp.expand(ser*(1-z)**(D+1)),z).all_coeffs()[::-1][:D+1]
      sgn=[sp.sign(x) for x in h if x!=0]
      ok=all(s==sgn[0] for s in sgn) if sgn else False
      res[o][ok]+=1
      if not ok: ex.setdefault(o,line.strip()[:120]+' h='+str(h))
for o in offs: print('c0=',o,dict(res[o]),ex.get(o,''))

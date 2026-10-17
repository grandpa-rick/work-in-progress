import sympy as sp, itertools, sys
from engine import subsets, geq, W, q
N=int(sys.argv[1]); w=sp.symbols('w1:%d'%(N+1))
S=subsets(N)
for M in S:
  for a,b in itertools.combinations(M,2):
    Mp=tuple(x for x in M if x not in (a,b))
    sub={w[b-1]:1/w[a-1]}
    rs=[]
    for I in S:
      if not geq(Mp,I): continue
      x=sp.cancel(W('C',N,I,M,q)[0].subs(sub)); y=sp.cancel(W('C',N,I,Mp,q)[0].subs(sub))
      rs.append((I,sp.factor(sp.cancel(x/y))))
    print(M,(a,b),Mp, 'ratios equal:',len(set(r for _,r in rs))==1, rs[0][1] if len(set(r for _,r in rs))==1 else rs, flush=True)

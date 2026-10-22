import sys, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di
def G(k,a):
    f={}
    for S in itertools.product((0,1),repeat=k):
        t=tuple(-a-s for s in S); f[t]=f.get(t,0)+1
    for i in range(k,0,-1): f=Di(f,i,k)
    return f
def Fb(k,b):
    f={}
    for S in itertools.product((0,1),repeat=k):
        t=tuple(-b-s for s in S); f[t]=f.get(t,0)+1
    return f
for k in range(1,7):
  for a in range(1,5):
    g=G(k,a)
    basis=[('F',b,Fb(k,b)) for b in range(-a-1,a+2)]+[('e',b,{tuple([-b]*k):1}) for b in range(-a-1,a+3)]
    keys=sorted(set(g)|set().union(*[set(B[2]) for B in basis]))
    M=sp.Matrix([[B[2].get(kk,0) for B in basis] for kk in keys]); v=sp.Matrix([g.get(kk,0) for kk in keys])
    try:
        sol,params=M.gauss_jordan_solve(v)
        sol=sol.subs({p:0 for p in params})
        print('k',k,'a',a,'IN SPAN:', {f'{B[0]}{B[1]}':c for B,c in zip(basis,sol) if c!=0})
    except ValueError:
        print('k',k,'a',a,'not in span')

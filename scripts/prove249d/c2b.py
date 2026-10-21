import sys, itertools, sympy as sp
from dem import chi as chiC
from demB import chiB
sys.path.insert(0,'../prove249c'); from eng import leq
n=int(sys.argv[1]); x=sp.Symbol('m'); c=sp.Symbol('c')
subs=[tuple(k) for r in range(n+1) for k in itertools.combinations(range(1,n+1),r)]
PC={};PB={}
for k in subs:
    d=sum(k)
    hc=sp.interpolate([(m,chiC(set(k),n,m)[0]) for m in range(d+2)],x)
    hb=sp.interpolate([(m,chiB(set(k),n,m)) for m in range(d+2)],x)
    PC[k]=sp.expand((-1)**d*hc.subs(x,-c/2-1)); PB[k]=sp.expand((-1)**d*hb.subs(x,-c-1))
D=n*(n+1)//2
for k in subs:
    below=[v for v in subs if leq(list(v),list(k))]
    cs=sp.symbols('b0:%d'%len(below))
    eq=sp.Poly(PC[k]-sum(ci*PB[v] for ci,v in zip(cs,below)),c).all_coeffs()
    sol=sp.linsolve(eq,cs)
    print(k,sol)

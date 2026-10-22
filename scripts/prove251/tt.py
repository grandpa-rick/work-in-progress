from pq import *
sys.path.insert(0,'.')
from gen import counter
from orderpoly import shifted, boxless
def T_of(k):
    P=shifted(sorted(k,reverse=True)); n=len(P); idx={b:i for i,b in enumerate(P)}
    F=counter(n, lambda a,b: boxless(P[a],P[b]), frozenset())
    return sp.expand(sp.interpolate([(N,F(N)) for N in range(n+2)],x))
m=int(sys.argv[1])
for r in range(1,m+1):
  for k in itertools.combinations(range(1,m+1),r):
    S,E,O,p,q=polys(k); T=T_of(k)
    print(k, 'T-E(N)-E(N-1)=', sp.factor(T-E-E.subs(x,x-1)), '   T-E-O(N-1)=',sp.factor(T-E-O.subs(x,x-1)))

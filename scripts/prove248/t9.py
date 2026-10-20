from psi import *
from t2 import lead0
import sys
def P2pred(s,m,Mm):
    tot=0
    for i in range(s-1):
        for T in itertools.combinations(range(1,s),i):
            mon=sp.Mul(*[w(k) for k in T])
            tot+=(-1)**i*mon*Psi((1,)+tuple(range(s,s+m-1)),(s-1-i,)+(s-1,)*(m-1),Mm)
    W=sp.Mul(*[w(k) for k in range(1,s)])
    tot+=(-1)**(s-1)*W*(Psi(tuple(range(s,s+m-1)),(s-1,)*(m-1),Mm) if m>1 else 1)
    return (-1)**m*tot
for s,m,M in [(2,2,(2,3,4)),(2,2,(1,2,3,5)),(3,2,(3,4,5)),(2,3,(2,3,4,5))]:
    n=M[-1]; u=w(n); Mm=M[:-1]
    o,c=lead0(Psi(tuple(range(s,s+m)),(s,)*m,M),u)
    print(s,m,M,'order',o,'(pred',-m,')',sp.simplify(c-P2pred(s,m,Mm)),flush=True)

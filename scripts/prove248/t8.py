# check P1' numerically: Psi_{(s),s}(M) * (-u) -> sum
from psi import *
from t2 import lead0
for s,M in [(2,(2,3,4)),(3,(3,4,5)),(2,(1,2,4))]:
    n=M[-1]; u=w(n); Mm=M[:-1]
    o,c=lead0(Psi((s,),(s,),M),u)
    pred=0
    for i in range(s):
        for T in itertools.combinations(range(1,s),i):
            mon=sp.Mul(*[w(k) for k in T])
            pred+= -(-1)**i*mon*(Psi((1,),(s-1-i,),Mm) if i<s-1 else 1)
    print(s,M,o,sp.simplify(c-pred))

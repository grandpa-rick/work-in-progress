# top u-degree coefficient of F_lam(M), n in M: compare with Pi_{M-} w^2 F_{lam-}(M-)
from fractions import Fraction as Fr
import itertools,random,sys
sys.path.insert(0,'../prove246'); from fast import S, Fv, subsets, geq
R=random.Random(1)
n=5; base=[Fr(R.randint(2,40),R.randint(2,40)) for _ in range(n-1)]
def F(lam,M,uval):
    w=[S.const(v) for v in base]+[S.const(uval)]
    return Fv(lam,M,w,1,NN=0).c[0]
for lam in [(1,),(3,),(4,),(1,2),(2,3),(3,4),(2,4),(1,2,3),(2,3,4),(1,3,4)]:
    rats=set()
    for M in subsets(n):
        if n not in M or not geq(M,lam): continue
        Mm=tuple(m for m in M if m!=n); lm=lam[:-1]
        U=Fr(10**12); v=F(lam,M,U)/U**lam[-1]
        ref=F(lm,Mm,Fr(1)) if geq(Mm,lm) else 0
        for m in Mm: ref*=base[m-1]**2
        rats.add(round(float(v/ref),6) if ref!=0 else ('ref0',round(float(v),6)))
    print(lam,rats)

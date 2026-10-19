from fractions import Fraction as Fr
import itertools,random,sys
sys.path.insert(0,'../prove246'); from fast import S, Fv, subsets, geq
R=random.Random(1)
n=5; base=[Fr(R.randint(2,40),R.randint(2,40)) for _ in range(n-1)]
def F(lam,M,uval):
    w=[S.const(v) for v in base]+[S.const(uval)]
    return Fv(lam,M,w,1,NN=0).c[0]
def Fm(lam,M):
    w=[S.const(v) for v in base]
    return Fv(lam,M,w,1,NN=0).c[0] if geq(M,lam) else 0
cands=[l for l in subsets(n-1)]
for lam in [(1,),(1,2),(2,3),(3,4),(1,3),(1,2,3),(2,3,4),(1,3,4),(1,2,4)]:
    U=Fr(10**15); psi={}
    for M in subsets(n):
        if n not in M or not geq(M,lam): continue
        Mm=tuple(m for m in M if m!=n); psi[Mm]=F(lam,M,U)/U**lam[-1]
    found=[]
    for e in range(-6,7):
        for lp in cands:
            rats=set()
            for Mm,v in psi.items():
                ref=Fm(lp,Mm)
                for m in Mm: ref*=base[m-1]**e
                if ref==0: rats.add('z' if abs(v)<1e-6 else 'bad')
                else: rats.add(round(float(v/ref),5))
            if 'bad' not in rats and len([r for r in rats if r!='z'])==1: found.append((e,lp,rats))
    print(lam,found[:4])

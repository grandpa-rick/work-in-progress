from fractions import Fraction as Fr
import random,sys
sys.path.insert(0,'../prove246'); from fast import S, Fv, subsets, geq
for seed in (1,2):
  R=random.Random(seed)
  n=5; base=[Fr(R.randint(2,40),R.randint(2,40)) for _ in range(n-1)]
  def F(lam,M,uval):
    w=[S.const(v) for v in base]+[S.const(uval)]
    return Fv(lam,M,w,1,NN=0).c[0]
  def Fm(lam,M):
    w=[S.const(v) for v in base]
    return Fv(lam,M,w,1,NN=0).c[0] if geq(M,lam) else 0
  def pw(M,e):
    r=Fr(1)
    for m in M: r*=base[m-1]**e
    return r
  for lam in [(1,2),(2,3),(3,4),(1,2,3),(2,3,4),(1,3,4)]:
    U=Fr(10**15); D={}
    for M in subsets(n):
        if n not in M or not geq(M,lam): continue
        Mm=tuple(m for m in M if m!=n)
        D[Mm]=F(lam,M,U)/U**lam[-1] - pw(Mm,2)*Fm(lam[:-1],Mm)
    found=[]
    for e in range(-6,7):
        for lp in subsets(n-1):
            rats=set()
            for Mm,v in D.items():
                ref=Fm(lp,Mm)*pw(Mm,e)
                if ref==0: rats.add('z' if abs(v)<1e-9 else 'bad')
                else: rats.add(round(float(v/ref),6))
            nz=[r for r in rats if r!='z']
            if 'bad' not in rats and len(nz)==1 and nz[0]!=0: found.append((e,lp,nz[0]))
    print(seed,lam,found[:5])

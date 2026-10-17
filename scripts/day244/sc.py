import sys; sys.path.insert(0,'../day242/q484')
from engine import *
from compare import koncki,t,s,y
w1,w2=sp.symbols('w1 w2'); d={w1:1/s,w2:1/t}
pts=['13','14','23','24']   # increasing BB = linear ext. for forward substitution (support: class ij nonzero at kl<=ij)
def basis(src,lamc):
  if src=="K": c,c2,al=lamc
  if src=='K': return {(cl,p):sp.sympify(koncki(cl,p,c,c2,al)) for cl in pts for p in pts}
  lab={():'24',(1,):'23',(2,):'14',(1,2):'13'}
  return {(lab[I],lab[M]):W('C',2,I,M)[0].subs(d).subs(q,-y) for I in subsets(2) for M in subsets(2)}
def sc(B):
  res={}
  for a in pts:
    for b in pts:
      c={}
      for p in pts[::-1]:  # top-down: f(p)=sum_{nu>=p} c_nu B[nu,p]
        rhs=B[a,p]*B[b,p]-sum(c[k]*B[k,p] for k in c)
        c[p]=sp.factor(sp.cancel(rhs/B[p,p]))
      res[a,b]=c
  return res
def poles(e):
  den=sp.fraction(sp.factor(sp.cancel(e)))[1]
  return [f for f,m in sp.factor_list(den)[1] if f.has(t) or f.has(s)]
for nm,B in [('Koncki lam=1/4',basis('K',(1,1,1))),('Motegi',basis('M',None))]:
  R=sc(B); allp=set()
  for k,c in R.items():
    for nu,e in c.items():
      for f in poles(e): allp.add(f)
  print(nm,'nonmonomial denominator factors in t,s:',allp)
  print('   c^{13}_{24,24} =',R['24','24']['13'])

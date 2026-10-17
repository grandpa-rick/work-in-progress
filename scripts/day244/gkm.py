import sys; sys.path.insert(0,'../day242/q484')
from engine import *
from compare import koncki, lab, S, t, s, y
w1,w2=sp.symbols('w1 w2'); d={w1:1/s,w2:1/t}
curves=[('24','14','t2'),('23','13','t2'),('24','23','s2'),('14','13','s2'),('24','13','st'),('23','14','t/s')]
subs={'t2':[{t:1},{t:-1}],'s2':[{s:1},{s:-1}],'st':[{t:1/s}],'t/s':[{t:s}]}
pts=['24','23','14','13']
def gkm(col):
  col={k:sp.sympify(v) for k,v in col.items()}
  return [(a,b,wt) for a,b,wt in curves for sb in subs[wt] if sp.simplify((col[a]-col[b]).subs(sb))!=0]
for lam in [sp.Rational(1,4),sp.Rational(3,4),0]:
  c=sp.ceiling(lam); c2=sp.ceiling(2*lam); al=2*c-c2
  print('lam',lam,'alpha',al)
  for cl in pts:
    print('  Koncki',cl,gkm({p:koncki(cl,p,c,c2,al) for p in pts}))
for I in S:
  print('  Motegi',lab[I],gkm({lab[M]:W('C',2,I,M)[0].subs(d).subs(q,-y) for M in S}))

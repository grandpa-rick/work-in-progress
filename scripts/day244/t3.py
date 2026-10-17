import sys; sys.path.insert(0,'../day242/q484')
from engine import *
from compare import koncki, lab, S, t, s, y
w1,w2=sp.symbols('w1 w2')
d={w1:1/s,w2:1/t}
c,c2,al=1,1,1
mot=sp.factor(W('C',2,(),(1,2))[0].subs(d).subs(q,-y))
kon=koncki('24','13',c,c2,al)
print('Motegi entry',mot); print('Koncki entry',sp.factor(kon))
print('Koncki-Motegi',sp.factor(sp.expand(kon-mot)))
r=sp.cancel(mot/kon); print('ratio on st=1:',sp.factor(r.subs(t,1/s)))
x=1/s; g=(1+y)**2/((1+y*x**2)*(1+y/x**2)); print('equals g(x=w1)?',sp.simplify(r.subs(t,1/s)-g)==0)
# GKM check of Koncki classes: curves and weights
pts=['24','23','14','13']
TW={'24':[t**2,s**2,t*s],'23':[t**2,1/s**2,t/s],'14':[1/t**2,s**2,s/t],'13':[1/t**2,1/s**2,1/(t*s)]}
curves=[]
for i,a in enumerate(pts):
  for b in pts[i+1:]:
    for wa in TW[a]:
      if 1/wa in TW[b] or sp.simplify(1/wa) in [sp.simplify(z) for z in TW[b]]: curves.append((a,b,wa))
print('curves',curves)
def gkm(col):
  bad=[]
  for a,b,wa in curves:
    diff=sp.together(col[a]-col[b]); num=sp.numer(diff)
    # divisible by (1-wa)?  substitute along wa=1
    if wa==t*s or wa==1/(t*s) or wa==t/s or wa==s/t: sub={t:(1/s if wa in (t*s,1/(t*s)) else s)}
    elif wa in (t**2,1/t**2): sub={t:1}
    else: sub={s:1}
    if sp.simplify(diff.subs(sub))!=0:
      # t^2=1 also t=-1
      bad.append((a,b,wa))
    if wa in (t**2,1/t**2) and sp.simplify(diff.subs({t:-1}))!=0: bad.append((a,b,wa,'-1'))
    if wa in (s**2,1/s**2) and sp.simplify(diff.subs({s:-1}))!=0: bad.append((a,b,wa,'-1'))
  return bad
for cl in pts:
  col={p:koncki(cl,p,c,c2,al) for p in pts}
  print('Koncki class',cl,'GKM fails:',gkm(col))
  colm={lab[M]:W('C',2,I,M)[0].subs(d).subs(q,-y) for I in S if lab[I]==cl for M in S}
  print('Motegi class',cl,'GKM fails:',gkm(colm))

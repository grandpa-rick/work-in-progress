import sys; sys.path.insert(0,'../day242/q484')
from engine import *
t,s,y=sp.symbols('t s y')
S=subsets(2)
lab={():'24',(1,):'23',(2,):'14',(1,2):'13'}
def koncki(cls,pt,c,c2,al):  # c=ceil(lam), c2=ceil(2lam)
  K={('13','13'):(1-s**2)*(1-t**2)*(1-s*t),
     ('14','14'):(1-t**2)*(1-t/s)*(1+y/s**2),
     ('14','13'):(1-t**2)*(1-s*t)*(1+y)*s**(2-2*c),
     ('23','23'):(1-s**2)*(1+y/t**2)*(1+y*s/t),
     ('23','14'):(1-t**2)*(1+y/s**2)*(1+y)*(s/t)**(c2-1),
     ('23','13'):(1-s*t)*(1+y)*t**(1-c2)*((s*(1+y)+t+s*t*y) if al==0 else (1+s*t*(1+y)+s**2*y)),
     ('24','24'):(1+y/t**2)*(1+y/s**2)*(1+y/(s*t)),
     ('24','23'):(1+y/t**2)*(1+y*s/t)*(1+y)*s**(2-2*c),
     ('24','14'):(1+y/s**2)*(1+y*t/s)*(1+y)*t**(2-2*c),
     ('24','13'):(1+y)*(s*t)**(1-c2)*((s*t*(1+y)-y+s**2*y+t**2*y+s**2*t**2*y**2) if al==0 else (s*t*(y+y**2)+1+s**2*y+t**2*y-s**2*t**2*y))}
  return K.get((cls,pt),0)
w1,w2=sp.symbols('w1 w2')
subsd=dict(swap={w1:1/t,w2:1/s},std={w1:1/s,w2:1/t})
for lam in [sp.Rational(-1,3),0,sp.Rational(1,4),sp.Rational(1,2),sp.Rational(2,3),1,sp.Rational(5,4),sp.Rational(3,2)]:
  c=sp.ceiling(lam); c2=sp.ceiling(2*lam); al=2*c-c2
  for nm,d in subsd.items():
    rows=[];ok=0
    for I in S:
      for M in S:
        v,_=W('C',2,I,M); v=sp.factor(v.subs(d).subs(q,-y))
        k=koncki(lab[I],lab[M],c,c2,al)
        if v==0 and k==0: ok+=1; continue
        r=sp.factor(sp.cancel(v/k)) if k!=0 else 'inf'
        if r==1: ok+=1
        else: rows.append((lab[I],lab[M],r))
    print('lam',lam,'alpha',al,nm,'match',ok,'/16')
    for rr in rows: print('   ',rr)

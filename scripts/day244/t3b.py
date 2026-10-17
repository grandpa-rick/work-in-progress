import sys; sys.path.insert(0,'../day242/q484')
from engine import *
from compare import koncki,t,s,y
w1,w2=sp.symbols('w1 w2'); d={w1:1/s,w2:1/t}
mot=sp.factor(W('C',2,(),(1,2))[0].subs(d).subs(q,-y)); kon=koncki('24','13',1,1,1)
print('Kon-Mot =',sp.factor(sp.expand(kon-mot)))
r=sp.factor(sp.cancel((mot/kon).subs(t,1/s))); x=1/s
g=(1+y)**2/((1+y*x**2)*(1+y/x**2)); print('ratio on st=1:',r,' == g:',sp.simplify(r-g)==0)
print('Motegi(e13)/Motegi(e24) on st=1:',sp.factor(sp.cancel((mot/W('C',2,(),())[0].subs(d).subs(q,-y)).subs(t,1/s))))
E=(1-s**2)*(1-t**2)*(1-s*t)
print('(Kon-Mot)/(point class) =',sp.factor(sp.cancel((kon-mot)/E)))
for al in [0,1]:
  for c2 in [2*1-al]:
    lrr=(1-s*t)*(1+y)*t**(1-c2)*((1+y*s/t)*(1-s**2)*t**(1-al)/(1-s/t)+(1+y)*(1-t**2)*s**(1-al)/(1-t/s))
    pr=koncki('23','13',1,c2,al)
    print('alpha',al,'LRR 23|e13 =',sp.factor(sp.cancel(lrr)),' printed =',sp.factor(pr),' equal:',sp.simplify(lrr-pr)==0)

from eng import *
import random, math
random.seed(1)
n=3
w={k:Fr(random.randint(2,30),random.randint(2,30)) for k in range(1,n+1)}
def D(w):
    p=Fr(1)
    for i in range(1,n+1):
        p*=(w[i]-1/w[i])
        for j in range(i+1,n+1): p*=(w[j]-w[i])*(w[i]*w[j]-1)
    return p
for lam,mu,kap in [((),(),(1,2,3)),((),(1,),(1,2,3)),((),(2,),(2,3)),((1,),(1,),(1,3))]:
  c=5
  P=D(w)*Cfun(list(lam),list(mu),list(kap),c,n,w)
  for i in range(1,n+1):
    w2=dict(w); w2[i]=1/w[i]
    r=D(w2)*Cfun(list(lam),list(mu),list(kap),c,n,w2)/P
    e=math.log(abs(r))/math.log(float(w[i]))
    print(lam,mu,kap,i,round(e,4),r<0)

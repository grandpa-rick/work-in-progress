from eng import *
import random
random.seed(1)
n=3
w={k:Fr(random.randint(2,30),random.randint(2,30)) for k in range(1,n+1)}
def D(w):
    p=Fr(1)
    for i in range(1,n+1):
        p*=(w[i]-1/w[i])
        for j in range(i+1,n+1): p*=(w[j]-w[i])*(w[i]*w[j]-1)
    return p
for lam,mu,kap in [((),(),(1,2,3)),((),(1,),(1,2,3)),((),(2,),(2,3)),((1,),(1,),(1,3)),((),(1,2),(1,2,3))]:
  for c in [5,6]:
    P=D(w)*Cfun(list(lam),list(mu),list(kap),c,n,w)
    out=[]
    for i in range(1,n+1):
        w2=dict(w); w2[i]=1/w[i]
        P2=D(w2)*Cfun(list(lam),list(mu),list(kap),c,n,w2)
        out.append(P2/P if P else None)
    print(lam,mu,kap,c,P!=0,[float(o) if o is not None else None for o in out])

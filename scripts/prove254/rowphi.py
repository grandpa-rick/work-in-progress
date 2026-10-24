# one-row reduction: Phi(u) = E0 - sum_{L>=1} (-1/2)^L G_L for a chain of n boxes, alphabet [M], bounds v_b < u_b (b=1..k)
import itertools, sys
from fractions import Fraction as Fr
def fills(n,M,u,L,even0):
    # count strictly increasing v_L..v_{n-1} in [1,M], v_b<u_b (b<=len(u)), v_0 even if L==0 and even0
    cnt=0
    for vs in itertools.combinations(range(1,M+1),n-L):
        ok=True
        for i,v in enumerate(vs):
            b=L+i
            if 1<=b<=len(u) and v>=u[b-1]: ok=False;break
            if b==0 and even0 and v%2: ok=False;break
        if ok: cnt+=1
    return cnt
def Phi(n,M,u):
    return fills(n,M,u,0,True)-sum(Fr(-1,2)**L*fills(n,M,u,L,False) for L in range(1,n+1))
neg=[]; tot=0; zero=0
for n in range(1,7):
  for M in range(n,n+5):
    for k in range(0,n):
      for u in itertools.combinations(range(2,M+1),k):
        if k and u[0]%2: continue
        p=Phi(n,M,u); tot+=1
        if p<0: neg.append((n,M,u,p))
        if p==0: zero+=1
print('tested',tot,'negative',len(neg),'zero',zero); print(neg[:10])

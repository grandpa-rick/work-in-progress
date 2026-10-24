from rowphi import *
import itertools
from fractions import Fraction as Fr
# general bounds: arbitrary weakly increasing? arbitrary vectors u_b in [1,M+1] for b=1..n-1
for Moff in (-1,0,1,2):
  neg=0; tot=0; zero=0; ex=[]
  for n in range(1,6):
    M=n+Moff
    if M<1: continue
    for u in itertools.product(range(1,M+2),repeat=n-1):
        p=Phi(n,M,u); tot+=1
        if p<0: neg+=1; ex.append((n,M,u,p)) if len(ex)<3 else None
        if p==0: zero+=1
  print('M-n=',Moff,'tested',tot,'neg',neg,'zero',zero,ex)

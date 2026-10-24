import sys
exec(open('tprime.py').read().split("ok=bad=0")[0])
from cnt import count
from fractions import Fraction as Fr
neg=0
for k in KS:
    if not k or k[0]%2==0: continue
    a=(k[0]+1)//2
    terms=[]
    for r in KS:
        if not boxes(r)<=boxes(k): continue
        t=Tp(boxes(k)-boxes(r))
        if t==0: continue
        Nv=count(r,2*a-1,True)
        if Nv: terms.append((r,str(t),Nv))
    tot=sum(Fr(t)*Nv for r,t,Nv in terms)
    if tot<0: neg+=1
    print(k,'a',a,'s=',tot,terms)
print('negative',neg)

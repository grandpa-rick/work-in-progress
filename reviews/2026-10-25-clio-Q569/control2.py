"""Negative control + paper's protrusion example (Rick 2026-10-25)."""
import itertools
from control import *
from tiles import labelled, normalise
n,k=5,2
P=parts_in(k,n-k); copro={nu:coproduct(nu,k) for nu in P}
Ph=parts_in(k,n-k-1)
# negative control: hexR in place of hex against d
ok=tot=0
for l,m,nu in itertools.product(Ph,Ph,Ph):
    if sum(nu)>sum(l)+sum(m): continue
    tot+=1; d=sgn(sum(l)+sum(m)-sum(nu))*copro[nu].get((l,m),0)
    ok+= count_puzzles(KTW+[NT['hexR']],n,l,m,nu,k)==d
print('NEG control: (KTW+hexR) vs d on hex domain:',ok,'/',tot)
# refined set from sec 8, macro strings
def macro(s):
    import re
    return [(int(d), {'r':0,'b':1,'bb':2}[t]) for t,d in re.findall(r'\\(bb|b|r)(\d)',s)]
REF=[NT['tri0U'],NT['tri0D']]+[normalise(labelled(macro(s))) for s in
 [r'\b0\b2\b4',r'\b0\b2\bb4',r'\bb1\b3\b5',r'\b0\r1\b3\r4',r'\b2\r3\b5\r0',r'\bb1\r2\b4\r5',r'\bb1\r2\bb4\r5',
  r'\r1\b2\r3\b4\r5\b0',r'\r1\b2\r3\bb4\r5\b0',r'\r1\b2\r4\b4\b0']]
lam,mu,nu=(2,1),(4,2),(4,3,1)
c=coproduct(nu,3); d=sgn(sum(lam)+sum(mu)-sum(nu))*c.get((lam,mu),0)
print('paper example d =',d)
print(' hex side 7:',count_puzzles(KTW+[NT['hex']],7,lam,mu,nu,3))
print(' hex side 8:',count_puzzles(KTW+[NT['hex']],8,lam,mu,nu,3))
print(' refined side 7:',count_puzzles(REF,7,lam,mu,nu,3))
# refined set census on 2x2 ambient at side n=k+m=4
k2,m2=2,2; Pr=parts_in(k2,m2); ok=tot=0
cp={nu:coproduct(nu,k2) for nu in Pr}
for l,m,nu in itertools.product(Pr,Pr,Pr):
    if sum(nu)>sum(l)+sum(m): continue
    tot+=1; dd=sgn(sum(l)+sum(m)-sum(nu))*cp[nu].get((l,m),0)
    ok+= count_puzzles(REF,k2+m2,l,m,nu,k2)==dd
print('refined set, ambient 2x2, side 4:',ok,'/',tot)

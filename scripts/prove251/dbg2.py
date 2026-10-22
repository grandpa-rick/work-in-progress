from fractions import Fraction as Fr
import itertools, pf
from dbg import Jz
lam=[3,2,1]; N=5
tot=Fr(0)
for q in range(4):
    for S in itertools.combinations([1,2,3],q):
        eps=1
        for r in [1,2,3]:
            if r not in S: eps*= -(-1)**lam[r-1]
        Jhat=(Jz(S)*2**sum(lam[s-1]+1 for s in S)) if S else 1
        tot+=(-1)**q*eps*Jhat
        print(S,'Jhat',Jhat)
print('decoupled sum normalized',tot)
print('Hhat_direct',pf.Hhat_direct(lam,N))

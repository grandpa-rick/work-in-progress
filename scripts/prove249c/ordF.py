from eng import F, leq
from fractions import Fraction as Fr
import itertools, math, sys
n=int(sys.argv[1]); a=[Fr(k*k+3,11)*(-1)**k for k in range(1,n+1)]
U=list(range(1,n+1)); subs=[list(S) for r in range(n+1) for S in itertools.combinations(U,r)]
def order(I,M):
    vs=[]
    for e in (Fr(1,10**12),Fr(1,10**13)):
        w={i:1+a[i-1]*e for i in range(1,n+1)}; vs.append(F(I,M,w))
    if vs[0]==0: return None
    return round(math.log10(abs(float(vs[0]/vs[1]))))
bad=0; tot=0; diag=0
for I in subs:
  for M in subs:
    if not leq(I,M): continue
    o=order(I,M); tot+=1
    if o is None: continue
    if o<sum(I): bad+=1; print('VIOL',I,M,o,sum(I))
    if I==M and o!=sum(I): diag+=1; print('DIAG',I,o)
print('n',n,'pairs',tot,'violations ord F_I(M) < |I|:',bad,'diag mismatches',diag)

import sys
import sys, itertools
sys.path.insert(0,".")
sys.argv=[sys.argv[0],sys.argv[1]]
import gate6 as G
from fractions import Fraction as Fr
N=int(sys.argv[1])
for L in range(1,N+1):
  for c in itertools.combinations(range(1,N+1),L):
    k=tuple(sorted(c,reverse=True)); rp=G.rp(k); M0=2*((rp+1)//2)+3
    s=''.join('0' if v==0 else ('+' if v>0 else '-') for v in [(-1)**sum(k)*G.H(k,Fr(-M,2)) for M in range(1,2*N+4,2)])
    print(k,"r'",rp,'M0',M0,s)

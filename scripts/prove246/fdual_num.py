import sys, itertools, random; sys.path.insert(0,'.')
import fast
from fast import *
from fractions import Fraction as Fr
fast.K=1
N=int(sys.argv[1]); R=random.Random(int(sys.argv[2]))
vals=[Fr(R.randint(2,60),R.randint(2,60)) for _ in range(N)]
w=[S.const(v) for v in vals]
Ss=sorted(subsets(N),key=lambda T:(len(T),sum(T)))
comp=lambda M: tuple(c for c in range(1,N+1) if c not in M)
def F(I,M,wt): return Fv(I,M,w,wt,NN=0).c[0] if geq(M,I) else Fr(0)
def cc(I):
    r=Fr(1)
    for i in I:
        r*=(-1)**i
        for k in range(1,i+1): r*=vals[k-1]
    return r
A={(M,I):F(I,M,1) for M in Ss for I in Ss}
P0={(M,I):F(I,M,0) for M in Ss for I in Ss}
# inverse by triangular solve per column e_L
bad=0;tot=0
for L in Ss:
    c={}
    for M in Ss:
        v=(Fr(1) if M==L else Fr(0))-sum(c[I]*A[M,I] for I in c if A[M,I]!=0)
        c[M]=v/A[M,M]
    for nu in Ss:
        pred=cc(L)/cc(nu)*P0[comp(L),comp(nu)]/(A[L,L]*P0[comp(L),comp(L)])
        tot+=1; bad+= (pred!=c[nu])
print('N',N,'mismatch',bad,'of',tot)

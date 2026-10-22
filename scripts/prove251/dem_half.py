import sys, itertools
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import Di, word
from pq import polys, x
import sympy as sp
def chiw(kap,n,f):
    ws=word(set(kap),n)
    for i in reversed(ws): f=Di(f,i,n)
    return sum(f.values())
n=int(sys.argv[1]); bad=0; tot=0
for r in range(1,n+1):
  for kap in itertools.combinations(range(1,n+1),r):
    S,E,O,p,q=polys(kap)
    for M in (1,3,5,7,9):
        a=(M-1)//2
        tot_=0
        for Ssub in itertools.product((0,1),repeat=n):
            tot_+=chiw(kap,n,{tuple(-a-s for s in Ssub):1})
        pred=sp.Rational(tot_,2**n); act=E.subs(x,-M)
        tot+=1; bad+= pred!=act
        if pred!=act: print('MISMATCH',kap,M,pred,act)
print('checked',tot,'bad',bad)

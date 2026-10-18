import sys, itertools, math, random; sys.path.insert(0,'.')
import fast
from fast import *
from fractions import Fraction as Fr
def pointx(N,a,b,seed,x):
    R=random.Random(seed)
    vals=[Fr(R.randint(2,97),R.randint(2,97))+R.randint(0,3) for _ in range(N)]
    vals[a-1]=Fr(x); w=[S.const(v) for v in vals]
    w[b-1]=S(0,[Fr(1)/Fr(x)]*fast.K); return w
class BX(Basis):
    def __init__(s,N,a,b,seed,wt,x,NN):
        s.N=N; s.w=pointx(N,a,b,seed,x); s.wt=wt; s.cache={}; s.a=a; s.b=b; s.NN=NN
    def A(s,M,I):
        if (M,I) not in s.cache: s.cache[M,I]=Fv(I,M,s.w,s.wt,NN=s.NN) if geq(M,I) else S.const(0)
        return s.cache[M,I]
def logabs(f): return math.log(abs(f.numerator))-math.log(abs(f.denominator)) if f!=0 else None
if __name__=='__main__':
    N=int(sys.argv[1]); a,b=int(sys.argv[2]),int(sys.argv[3]); l=tuple(map(int,sys.argv[4].replace('_',''))); m=tuple(map(int,sys.argv[5].replace('_',''))); nu=tuple(map(int,sys.argv[6]))
    for wt in (0,1):
        vals=[]
        for k in (4,8,12):
            x=10**k
            B=BX(N,a,b,3,wt,x,N if wt else 0)
            r=B.resid(l,m,nu); vals.append((k,r))
        print('wt',wt,[ (k, float(r)) for k,r in vals])
        L=[logabs(r) for k,r in vals]
        if None not in L: print('   slopes in log10 x:',(L[1]-L[0])/(4*math.log(10)),(L[2]-L[1])/(4*math.log(10)))

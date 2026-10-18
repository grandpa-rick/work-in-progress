"""R(c) at w_C = 1 + s r_c, s->0 (x fixed). Is the limit finite, and polynomial in c?"""
import sys, itertools, random; sys.path.insert(0,'.')
import fast
from fast import *
from fractions import Fraction as Fr
class BW(Basis):
    def __init__(s,N,a,b,wt,x,sval,seed=5):
        R=random.Random(seed); s.N=N; s.wt=wt; s.cache={}; s.a=a; s.b=b
        vals=[1+sval*Fr(R.randint(1,50),R.randint(1,50)) for _ in range(N)]
        vals[a-1]=Fr(x); s.w=[S.const(v) for v in vals]; s.w[b-1]=S(0,[Fr(1)/Fr(x)]*fast.K)
    def A(s,M,I):
        if (M,I) not in s.cache: s.cache[M,I]=Fv(I,M,s.w,s.wt,NN=0) if geq(M,I) else S.const(0)
        return s.cache[M,I]
    def residc(s,lam,mu,nu,c):
        a,b=s.a,s.b; f={}
        for M in subsets(s.N):
            if a in M and b in M and geq(nu,M):
                wgt=S.const(1)
                for m in M: wgt=wgt*(s.w[m-1]**c)
                f[M]=s.A(M,lam)*s.A(M,mu)*wgt
        return s.coords(f,nu).coef(-1)
if __name__=='__main__':
    N=int(sys.argv[1]); a,b=int(sys.argv[2]),int(sys.argv[3]); l=tuple(map(int,sys.argv[4].replace('_',''))); m=tuple(map(int,sys.argv[5].replace('_',''))); nu=tuple(map(int,sys.argv[6])); wt=int(sys.argv[7])
    x=Fr(3)
    for c in range(-10,11,2):
        r=[float(BW(N,a,b,wt,x,Fr(1,10**k)).residc(l,m,nu,c)) for k in (6,9)]
        print(c, r)

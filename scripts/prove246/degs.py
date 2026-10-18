"""For a case (N,a,b,lam,mu,nu,wt), list E and, for each spectator m and direction (inf/0), the w_m-degree of
 A_{M''} := Res l_nu(e_{M''+ab}) * F_lam(M'')F_mu(M'')   (no weight)."""
import sys, itertools, random, math; sys.path.insert(0,'.')
import fast
from fast import *
from fractions import Fraction as Fr
class BV(Basis):
    def __init__(s,N,a,b,wt,vals):
        s.N=N; s.wt=wt; s.cache={}; s.a=a; s.b=b
        s.w=[S.const(v) for v in vals]; s.w[b-1]=S(0,[Fr(1)/vals[a-1]]*fast.K)
    def A(s,M,I):
        if (M,I) not in s.cache: s.cache[M,I]=Fv(I,M,s.w,s.wt,NN=0) if geq(M,I) else S.const(0)
        return s.cache[M,I]
    def terms(s,lam,mu,nu):
        a,b=s.a,s.b; out={}
        for M in subsets(s.N):
            if a in M and b in M and geq(nu,M):
                Mpp=tuple(c for c in M if c not in (a,b))
                P=(s.A(Mpp,lam)*s.A(Mpp,mu)).coef(0)
                if P==0: continue
                rho=s.coords({M:S.const(1)},nu).coef(-1)
                out[Mpp]=rho*P
        return out
def la(f): return math.log10(abs(f.numerator))-math.log10(abs(f.denominator))
def run(N,a,b,lam,mu,nu,wt,seed=7):
    R=random.Random(seed); base=[Fr(R.randint(2,97),R.randint(2,97))+R.randint(0,3) for _ in range(N)]
    res={}
    for m in range(1,N+1):
        if m in (a,b): continue
        for dirn,ks in (('inf',(6,12)),('0',(-6,-12))):
            T=[]
            for k in ks:
                v=list(base); v[m-1]=v[m-1]*Fr(10)**k
                T.append(BV(N,a,b,wt,v).terms(lam,mu,nu))
            res[m,dirn]={M:round((la(T[1][M])-la(T[0][M]))/(ks[1]-ks[0]),2) for M in T[0] if T[0][M]!=0 and T[1].get(M,0)!=0}
    return res
if __name__=='__main__':
    N=int(sys.argv[1]); a,b=int(sys.argv[2]),int(sys.argv[3]); l=tuple(map(int,sys.argv[4].replace('_',''))); m=tuple(map(int,sys.argv[5].replace('_',''))); nu=tuple(map(int,sys.argv[6])); wt=int(sys.argv[7])
    for k,v in run(N,a,b,l,m,nu,wt).items(): print(k,v)

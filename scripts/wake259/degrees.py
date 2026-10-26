# Wake 259 Task 1: verify rem:degrees  deg X_k = 2^{|k|-l(k)} e(P_k) = 2^{|k|-l(k)} deg X^OG_k.
# Three independent routes for deg X_k (LG(m,2m)):
#  (A) leading coeff of H_k(j)=chi(X_k,O(j)) from Demazure operators (prove249d/dem.py), times |k|!
#  (B) Chevalley formula: sum over saturated Bruhat chains in W/W_P from the point to X_k of prod <w.omega_m, beta^vee>
#  (C) 2^{|k|-l} * e(P_k), e by brute-force linear extensions of the shifted diagram; cross-checked with
#      the shifted hook-length formula, and (D) OG degree via type-B Chevalley chains (spin weight).
import sys, itertools, math
from fractions import Fraction as Fr
from functools import lru_cache
import sympy as sp
sys.path.insert(0,'../prove249d'); from dem import chi, word

def shifted(k):
    return [(a,b) for a in range(len(k)) for b in range(a,a+k[a])]
def linext(P):
    P=list(P); n=len(P)
    less={x:[y for y in P if y!=x and y[0]<=x[0] and y[1]<=x[1]] for x in P}
    @lru_cache(None)
    def f(S):
        if len(S)==n: return 1
        return sum(f(S|frozenset([x])) for x in P if x not in S and all(y in S for y in less[x]))
    return f(frozenset())
def shifted_hook(k):
    D=shifted(k); s=set(D); n=len(D); prod=1
    for (a,b) in D:
        arm=sum(1 for c in range(b+1,a+k[a]) )
        leg=sum(1 for r in range(a+1,len(k)) if (r,b) in s)
        # shifted hook: also continue into row b+1 if box (a,b) lies left of row b+1's start column... standard:
        extra = k[b+1] if b+1 < len(k) else 0
        prod*= arm+leg+1+extra
    return math.factorial(n)//prod

def chev(m, typ):
    # weights v in {+-1}^m (orbit of omega_m; for B use spin weight (1/2,...,1/2) -> v/2)
    if typ=='C':
        R=[]
        for i in range(m):
            for k in range(i+1,m):
                for s in (1,-1):
                    v=[0]*m; v[i]=1; v[k]=-s; R.append((tuple(v),tuple(v)))
            v=[0]*m; v[i]=2; c=[0]*m; c[i]=1; R.append((tuple(v),tuple(c)))
        sc=Fr(1)
    else:  # B_m, short roots e_i coroot 2e_i, weight = v/2
        R=[]
        for i in range(m):
            for k in range(i+1,m):
                for s in (1,-1):
                    v=[0]*m; v[i]=1; v[k]=-s; R.append((tuple(v),tuple(v)))
            v=[0]*m; v[i]=1; c=[0]*m; c[i]=2; R.append((tuple(v),tuple(c)))
        sc=Fr(1,2)
    pr=lambda a,b: sum(Fr(x)*y for x,y in zip(a,b))
    V=list(itertools.product([1,-1],repeat=m))
    dim=lambda v: sum(1 for (a,c) in R if pr(v,c)<0)
    @lru_cache(None)
    def deg(v):
        if dim(v)==0: return 1
        tot=0
        for (a,c) in R:
            p=pr(v,c)*sc
            # u = s_beta v, need u covered by v (dim one less); bond = |<u*sc,beta^vee>|
            u=tuple(int(x-pr(v,c)*y) if typ=='C' else int(x-pr(v,c)*y) for x,y in zip(v,a))
            if dim(u)==dim(v)-1: tot+=abs(pr(u,c)*sc)*deg(u)
        return tot
    return deg
if __name__=='__main__':
    M=int(sys.argv[1]); AMAX=int(sys.argv[2]); J=sp.Symbol('j')
    print('kappa | dim | (A) Demazure lead*|k|! | (B) Chevalley LG | (C) 2^{|k|-l} e(P) | e(P) | shifted-hook | (D) Chevalley OG | verdict')
    allok=True
    for m in range(1,M+1):
        dC=chev(m,'C'); dB=chev(m,'B')
        for L in range(m+1):
            for c in itertools.combinations(range(1,m+1),L):
                k=tuple(sorted(c,reverse=True)); n=sum(k)
                if n>AMAX: vals=None
                else: vals=[chi(set(k),m,t)[0] for t in range(n+2)]
                if vals is not None: assert chi(set(k),m,0)[1]==n
                P=None if vals is None else sp.Poly(sp.interpolate(list(zip(range(n+2),vals)),J),J) if n>0 else sp.Poly(1,J)
                A=None if P is None else P.LC()*math.factorial(n)
                kv=tuple(-1 if (m-i) in c else 1 for i in range(m))
                B=dC(kv); Dg=dB(kv)
                e=linext(shifted(k)); C=2**(n-len(k))*e; h=shifted_hook(k) if k else 1
                ok = (A is None or A==B) and B==C and e==h==Dg
                allok&=ok
                print(f'm={m} {k} | {n} | {A if P is not None else "skip"} | {B} | {C} | {e} | {h} | {Dg} | {"ok" if ok else "FAIL"}')
    # intrinsic convention check: X_(2,1) in LG(3,6) must be LG(2,4) (word avoids s_1)
    print('word for (2,1) in m=3 (simple reflections used):',sorted(set(word({2,1},3))),' word for (3):',sorted(set(word({3},3))))
    print('ALL PASS' if allok else 'SOME FAIL')

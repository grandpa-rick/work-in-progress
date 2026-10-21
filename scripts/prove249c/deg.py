# xi-degree of T_L / w_L^c via exact evaluation at w_i = a_i * t^xi_i, t = 2^K
import itertools, random, math
from fractions import Fraction as Fr
from eng import Psi, F, ell, leq
def deg_eval(fn, xi, a, K1=60, K2=120):
    vals=[]
    for K in (K1,K2):
        t=Fr(2)**K
        w={i+1: a[i]*t**xi[i] for i in range(len(xi))}
        v=fn(w)
        if v==0: return None
        vals.append(v)
    l1=math.log2(abs(vals[0].numerator))-math.log2(vals[0].denominator)
    l2=math.log2(abs(vals[1].numerator))-math.log2(vals[1].denominator)
    d=(l2-l1)/(K2-K1)
    assert abs(d-round(d))<1e-3, d
    return round(d)
def hL(lam,mu,kap,M,n):
    return lambda w: ell(kap,M,n,w)*F(lam,M,w)*F(mu,M,w)

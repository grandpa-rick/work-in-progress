# Phi_a(p_nu; x,y) = <T_a p_nu, p_x p_y> via adjoint formula + exact [z^mu]K. grade: computed
from nonsym import *
import itertools, functools
def I_sym(a, nu, negs):
    """CT[ p_nu(z) * Z * prod_{m in negs} p_m(1/z) * K ] exact."""
    tot = 0
    for pos in itertools.product(range(a), repeat=len(nu)+len(negs)):
        gam = [1]*a
        for r, p in zip(list(nu)+[-m for m in negs], pos): gam[p] += r
        if gam[0] > 0: continue
        tot += ctmono(tuple(gam))
    return sp.expand(tot)
def phia(a): return sp.prod([1-t**i for i in range(1, a+1)])
@functools.lru_cache(None)
def Phi(a, nu, x, y):
    return sp.factor(sp.cancel((1-t**x)*(1-t**y)/phia(a)*I_sym(a, nu, (x, y))))
if __name__ == '__main__':
    from fractions import Fraction as Fr
    from twopoint_helpers import phi_engine
    for a, nu, x, y in [(2,(2,),2,2),(2,(1,1),3,1),(3,(2,1),4,2),(3,(3,2),5,3),(2,(3,),3,2)]:
        v = Phi(a, nu, x, y)
        print(a, nu, x, y, v, float(v.subs(t, sp.Rational(3,10))), float(phi_engine(a, nu, x, y, Fr(3,10))))

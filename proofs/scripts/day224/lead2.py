from ops import *
import sympy
ts = sympy.symbols('t')
def emul(f, g):  # product in e-basis dicts (keys = partitions)
    return mul(f, g)
def Dder(a, F, t):  # derivation D_a on an e-basis dict, via Thm 7.1
    out = {}
    for key, v in F.items():
        for i, x in enumerate(key):
            rest = key[:i] + key[i+1:]
            for k2, w in Mthm(a, x, t).items():
                kk = tuple(sorted(rest + k2, reverse=True)); out[kk] = out.get(kk, 0) + v*w
    return {k: v for k, v in out.items() if v}
def second(a, b, c, t):
    """[(s-1)^2] E_aE_bE_c(1) in e-basis = E2_a(e_b e_c) + D_a D_b e_c + e_a E2_b(e_c)."""
    X = toe(E2prod(a, b, c, t), t)
    Y = Dder(a, Mthm(b, c, t), t)
    Z = mul({(a,): Fr(1)}, toe(E2op(b, c, t), t))
    return add(add(X, Y), Z)
def pieces(a, b, c, t):
    return dict(Gam=toe(Gam(a, b, c, t), t), Ebc=mul({(b,): Fr(1)}, toe(E2op(a, c, t), t)), Ecb=mul({(c,): Fr(1)}, toe(E2op(a, b, t), t)),
                DD=Dder(a, Mthm(b, c, t), t), eE=mul({(a,): Fr(1)}, toe(E2op(b, c, t), t)))
def interp(fn, npts=40):
    """fn(t)->Fraction, polynomial in t; interpolate."""
    xs = list(range(2, 2+npts)); ys = [fn(Fr(x)) for x in xs]
    p = sympy.interpolate(list(zip(xs, ys)), ts)
    assert sympy.Poly(p, ts).degree() < npts - 5, 'degree too high'
    return sympy.factor(p)

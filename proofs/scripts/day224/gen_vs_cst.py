from second_gen import *
import pickle, sympy
C = {}
for n in range(2, 7):
    f = '../wake222/gprime/cst_n%d.pkl' % n if n < 6 else '../wake222/gprime/cst_n6_skip1n.pkl'
    C.update(pickle.load(open(f, 'rb')))
s, ts, e = sympy.symbols('s t e')
ok = bad = 0
for lam in set(k[0] for k in C):
    if len(lam) < 4: continue
    S2 = second_general(lam, Fr(3))[2]
    for (l2, mu), c in C.items():
        if l2 != lam: continue
        v = sympy.expand(sympy.sympify(c).subs({ts: 3, s: 1+e})).coeff(e, 2)
        g = sympy.Rational(S2.get(mu, 0)) == v; ok += g; bad += not g
        if not g: print('FAIL', lam, mu, v, S2.get(mu, 0))
print(ok, bad)

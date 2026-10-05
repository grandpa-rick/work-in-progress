from closed import *
import pickle
C = {}
for n in range(2, 7):
    f = '../wake222/gprime/cst_n%d.pkl' % n if n < 6 else '../wake222/gprime/cst_n6_skip1n.pkl'
    C.update(pickle.load(open(f, 'rb')))
ok = bad = 0
for (lam, mu), c in C.items():   # star2 check
    if len(lam) == 2:
        v = star2(*lam).get(mu, 0)
        g = sympy.simplify(sympy.sympify(c) - v) == 0; ok += g; bad += not g
        if not g: print('star2 FAIL', lam, mu, c, sympy.factor(v))
print('star2 vs engine', ok, bad)
ok = bad = 0
for (lam, mu), c in C.items():
    if len(lam) == 3 and mu == (sum(lam)-1, 1):
        for perm in set(__import__('itertools').permutations(lam)):
            v = cN1(perm)
            g = sympy.simplify(sympy.sympify(c) - v) == 0; ok += g; bad += not g
            if not g: print('cN1 FAIL', perm, mu, sympy.factor(c), sympy.factor(v))
print('cN1 vs engine', ok, bad)

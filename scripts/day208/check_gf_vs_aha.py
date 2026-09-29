"""Check symbolic Gamma_k (gf_recursion.py) against direct AHA at numeric points, exact rationals.
LHS: sum_{a,b<=m} z^a w^b * [t^{-C(k,2)} e_k(Y).(e_a e_b)](X)   (engine: ek_two_col.ekY, flint, m variables)
RHS: Gamma_k(z,w) evaluated at X (e_n(X), E(t^i z)=prod(1+t^i z X_l))."""
import sys, pickle, random
from fractions import Fraction as Fr
import sympy as sp
from ek_two_col import ekY
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA
from gf_recursion import s, t, z, w
K = int(sys.argv[1]); res = pickle.load(open(sys.argv[2], 'rb'))
random.seed(5)
def evalp(A, P, X, sv, tv):
    tot = Fr(0)
    for mon, cc in P.to_dict().items():
        v = Fr(int(cc))
        for i in range(A.m): v *= X[i]**int(mon[i])
        tot += v*sv**int(mon[A.m])*tv**int(mon[A.m+1])
    return tot
def ev(X, r):
    if r < 0 or r > len(X): return Fr(0)
    import itertools, math
    return sum((math.prod(c) for c in itertools.combinations(X, r)), Fr(0)) if r else Fr(1)
ok = True
for m in range(1, 5):
    A = AHA(m)
    for k in range(1, K+1):
        ops = {(a, b): ekY(A, A.e(a)*A.e(b), k) if k <= m else None for a in range(m+1) for b in range(m+1)}
        for trial in range(2):
            X = [Fr(random.randint(1, 30), random.randint(31, 60)) for _ in range(m)]
            sv, tv = Fr(random.randint(2, 9), 11), Fr(random.randint(2, 9), 13)
            zv, wv = Fr(random.randint(1, 20), 7), Fr(random.randint(21, 40), 9)
            lhs = sum((zv**a*wv**b*evalp(A, P, X, sv, tv) for (a, b), P in ops.items() if P is not None), Fr(0))
            rhs = Fr(0)
            import math
            for (lam, i, j), R in res[k].items():
                Rv = sp.Rational(R.subs({s: sp.Rational(sv.numerator, sv.denominator), t: sp.Rational(tv.numerator, tv.denominator),
                                        z: sp.Rational(zv.numerator, zv.denominator), w: sp.Rational(wv.numerator, wv.denominator)}))
                Rv = Fr(int(Rv.p), int(Rv.q))
                el = math.prod((ev(X, n) for n in lam), start=Fr(1))
                Ez = math.prod((1 + tv**i*zv*x for x in X), start=Fr(1)); Ew = math.prod((1 + tv**j*wv*x for x in X), start=Fr(1))
                rhs += Rv*el*Ez*Ew
            good = (lhs == rhs); ok &= good
            print(f'm={m} k={k} trial={trial}: {"OK" if good else "FAIL"}', flush=True)
print('ALL OK' if ok else 'SOME FAIL')

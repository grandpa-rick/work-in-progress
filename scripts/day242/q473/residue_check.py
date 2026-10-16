"""Day 242 Q473: independent exact check of Clio Prop 18 (iterated residues at x_t=(q y_{kappa_t})^{-1}).
Exact SymPy over Q at a random rational parameter point. G built from printed (5.1), c->oo.
Residue taken as lim_{x->p} (x-p)*F on the WHOLE G (sum of terms, cancelled), not term-by-term,
so no use is made of Clio's 'only phi has the pole' argument."""
import sympy as sp, random
from itertools import combinations, permutations
from gdmw_modp import subsets, preceq

R = random.Random(9)
fr = lambda: sp.Rational(R.randint(2, 97), R.randint(101, 199))
q, a = fr(), fr()
y = {j: fr() * R.choice([1, 3, 5, 7]) for j in range(1, 6)}
zeta = lambda x: x * (1 - a * x) / (x - a)
hinf = lambda x: a * (1 - x**2) / (a - x)


def phi(m, x):
    v = (1 - q) * x * y[m] / (1 - q * x * y[m])
    for j in range(1, m):
        v *= (1 - x * y[j]) / (1 - q * x * y[j])
    return v


def G(kappa, xs):
    d, L = len(kappa), len(xs)
    tot = 0
    for K in combinations(range(L), d):
        Kb = [j for j in range(L) if j not in K]
        pre = sp.Integer(1)
        for j in Kb:
            pre *= zeta(xs[j])
        for i in K:
            pre *= hinf(xs[i])
            for j in Kb:
                pre *= (xs[j] - q * xs[i]) / (xs[j] - xs[i]) * (1 - xs[i] * xs[j]) / (1 - q * xs[i] * xs[j])
        for ii in range(d):
            for jj in range(ii + 1, d):
                pre *= (1 - xs[K[ii]] * xs[K[jj]]) / (1 - q * xs[K[ii]] * xs[K[jj]])
        s = 0
        for sg in permutations(range(d)):
            t = sp.Integer(1)
            for i in range(d):
                for j in range(i + 1, d):
                    xa, xb = xs[K[sg[j]]], xs[K[sg[i]]]
                    t *= (xa - q * xb) / (xa - xb)
            for i in range(d):
                t *= phi(kappa[i], xs[K[sg[i]]])
            s += t
        tot += pre * s
    return tot


def Rk(k, kp, L):
    d = len(k)
    X = sp.symbols('x1:%d' % (d + 1))
    free = [fr() for _ in range(L - d)]
    F = sp.cancel(sp.together(G(kp, list(X) + free)))
    for t in range(d):
        p = 1 / (q * y[k[t]])
        F = sp.cancel((X[t] - p) * F)
        F = sp.cancel(F.subs(X[t], p))
    return F


for L, maxel in ((2, 4), (3, 4)):
    S = subsets(maxel, L)
    mis = z = nz = dg = 0
    rows = [k for k in S if k]
    for k in rows:
        for kp in S:
            v = Rk(k, kp, L)
            pred = preceq(k, kp)  # Clio: R_k(G_kp) != 0 iff kp >= k
            mis += (v != 0) != pred
            z += v == 0
            nz += v != 0
            if kp == k:
                dg += v != 0
    print(f"L={L} max={maxel}: cells {len(rows)*len(S)}, zero {z}, nonzero {nz}, "
          f"mismatch vs 'nonzero iff kp >= k': {mis}, diagonal nonzero {dg}/{len(rows)}", flush=True)

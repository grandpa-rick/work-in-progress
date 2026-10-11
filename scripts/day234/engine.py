# Day 234 ground-truth engine (written fresh this session; does NOT import any earlier engine).
# Implements ONLY Hikita's subset formula (FPSAC Prop 2.1, as printed):
#   e_k * F = E_k F = sum_{|A|=k} c_A X_A F(x_{A^c}, s x_A),  c_A = prod_{i in A, j notin A} (x_i - t x_j)/(x_i - x_j)
# in N variables, with symmetric functions stored in the e-basis {partition tuple: poly in Q[s,t]}.
# Evaluation strategy: numeric integer points x in Z^N, s,t symbolic; c_A is then a polynomial in t.
# E_k G is recovered in the e-basis by solving the (numeric) evaluation system at p(D) points.
import itertools, random
from fractions import Fraction
from functools import lru_cache
from sympy.polys.rings import ring
from sympy import QQ, Matrix, Rational

R, s, t = ring("s,t", QQ)

def partitions(n, maxpart=None, maxlen=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    if maxlen == 0: return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k, None if maxlen is None else maxlen - 1):
            yield (k,) + rest

def esym(vals, k):
    """elementary symmetric e_0..e_k of a list of ring elements / numbers"""
    e = [R(1)] + [R(0)] * k
    for v in vals:
        for j in range(k, 0, -1):
            e[j] = e[j] + e[j - 1] * v
    return e

class Engine:
    def __init__(self, N, seed=1234):
        self.N = N
        self.rng = random.Random(seed)
        self._pts = {}

    def points(self, D):
        if D in self._pts: return self._pts[D]
        basis = list(partitions(D, maxlen=self.N))
        while True:
            pts = [[self.rng.randint(-40, 40) for _ in range(self.N)] for _ in basis]
            ok = all(len(set(p)) == self.N and 0 not in p for p in pts)
            if not ok: continue
            M = Matrix([[self.eval_e_num(mu, p) for mu in basis] for p in pts])
            if M.det() != 0:
                Minv = M.inv()
                self._pts[D] = (basis, pts, Minv)
                return self._pts[D]

    @staticmethod
    def eval_e_num(mu, x):
        e = [1] + [0] * (max(mu) if mu else 0)
        for v in x:
            for j in range(len(e) - 1, 0, -1):
                e[j] += e[j - 1] * v
        r = 1
        for m in mu: r *= e[m]
        return r

    def eval_G(self, G, vals):
        """G in e-basis, vals = list of ring elements (the substituted variables)"""
        if not G: return R(0)
        mx = max((max(mu) if mu else 0) for mu in G)
        e = esym(vals, mx)
        tot = R(0)
        for mu, c in G.items():
            term = c
            for m in mu: term = term * e[m]
            tot += term
        return tot

    def cA(self, A, x):
        num = R(1); den = Fraction(1)
        for i in A:
            for j in range(self.N):
                if j in A: continue
                num = num * (x[i] - t * x[j])
                den *= (x[i] - x[j])
        return num * R(QQ(den.denominator) / QQ(den.numerator))

    def solve(self, D, values):
        basis, pts, Minv = self.points(D)
        out = {}
        for a, mu in enumerate(basis):
            acc = R(0)
            for i in range(len(pts)):
                m = Minv[a, i]
                if m != 0:
                    acc += values[i] * R(QQ(int(m.p), int(m.q)))
            if acc != 0: out[mu] = acc
        return out

    def E(self, k, G, deg):
        """E_k G, G homogeneous of degree deg in e-basis"""
        D = k + deg
        basis, pts, _ = self.points(D)
        vals = []
        for x in pts:
            tot = R(0)
            for A in itertools.combinations(range(self.N), k):
                XA = 1
                for i in A: XA *= x[i]
                sub = [R(x[j]) if j not in A else s * x[j] for j in range(self.N)]
                tot += self.cA(A, x) * XA * self.eval_G(G, sub)
            vals.append(tot)
        return self.solve(D, vals)

    def T(self, k, f, deg):
        """T_k f := sum_{|A|=k} c_A X_A f(x_A), f in e-basis (as symmetric fn in k vars)"""
        D = k + deg
        basis, pts, _ = self.points(D)
        vals = []
        for x in pts:
            tot = R(0)
            for A in itertools.combinations(range(self.N), k):
                XA = 1
                for i in A: XA *= x[i]
                sub = [R(x[i]) for i in A]
                tot += self.cA(A, x) * XA * self.eval_G(f, sub)
            vals.append(tot)
        return self.solve(D, vals)

    def estar(self, lam):
        """e*_lambda = E_{l1} ... E_{l_ell}(1)"""
        G = {(): R(1)}; deg = 0
        for k in reversed(lam):
            if k == 0: continue
            G = self.E(k, G, deg); deg += k
        return G

def mul_e(G, H):
    out = {}
    for a, c in G.items():
        for b, d in H.items():
            mu = tuple(sorted(a + b, reverse=True))
            out[mu] = out.get(mu, R(0)) + c * d
    return {k: v for k, v in out.items() if v != 0}

def sub_s(poly, sval):
    """substitute s -> sval (ring element or number)"""
    return poly.compose(s, sval) if not isinstance(sval, int) else poly.compose(s, R(sval))

def coeff_sm1(poly, p):
    """[(s-1)^p] poly"""
    q = poly.compose(s, s + 1)
    out = R(0)
    for (i, j), c in q.terms():
        if i == p: out += c * t**j
    return out

def val_sm1(poly):
    if poly == 0: return None
    q = poly.compose(s, s + 1)
    return min(i for (i, j), c in q.terms())

def val_s(poly):
    if poly == 0: return None
    return min(i for (i, j), c in poly.terms())

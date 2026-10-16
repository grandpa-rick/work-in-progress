"""Day 242: check Clio's suggested corollary (review 2026-10-15, Finding 3).

(i)  F_n(t^{-n}) = 0 for all n >= 1   (linear factor (1 - t^n w) of F_n; m = 0 case of Thm 1)
(ii) [z^{-N}] T_k = sum_{b=0}^{k-N} s^b F_{k-b}(t^{-N-b}) e_b e_{k-N-b} = 0,  1 <= N <= k,
     e_1, e_2, ... independent indeterminates (identity in Lambda).

Definitions transcribed from longversion.tex Sec. 3 (c(n,j), alpha_j, F_n).
Arms:
  A. (i) exact in QQ(s,t), n = 1..12, with control F_n(t^{-(n+1)}) != 0.
  B. (ii) exact in QQ(s,t), k = 1..8, every N in [1,k] (generic interior N included),
     computed by TWO mechanisms: the closed coefficient formula, and direct Laurent
     expansion of the generating function sum_{j+b<=k} s^b t^{-kj} c(k-b,j) e_b z^{b-k} E(t^j z).
     Positive control: [z^0]T_k = e_k (= e_k * 1).  Negative control: [z^{+1}]T_k != 0.
  C. Lesion: perturb c(n,j) (s t^{n-j} -> s t^{n-j+1}); (i) and (ii) must FAIL.
  D. Transcription sanity: Pieri rule vs. the subset operator E_k (Thm subset) in m variables,
     exact rationals, m <= 4, k <= 3, r <= 3.
"""
import itertools, sys
from collections import defaultdict
import sympy as sp
from sympy import QQ
from sympy.polys.rings import ring

R, sR, tR = ring('s,t', QQ)
sS, tS = sp.symbols('s t')

class Fr:
    """Unreduced fraction num/den over QQ[s,t] (no gcds: avoids sympy heugcd failures).
    Zero test = numerator is the zero polynomial, which is exact."""
    __slots__ = ('n', 'd')
    def __init__(self, n, d=None):
        self.n = n if not isinstance(n, int) else R(n)
        self.d = R(1) if d is None else d
    @staticmethod
    def c(x):
        return x if isinstance(x, Fr) else Fr(R(x))
    def __add__(a, b):
        b = Fr.c(b)
        if a.d == b.d: return Fr(a.n + b.n, a.d)
        return Fr(a.n * b.d + b.n * a.d, a.d * b.d)
    __radd__ = __add__
    def __neg__(a): return Fr(-a.n, a.d)
    def __sub__(a, b): return a + (-Fr.c(b))
    def __rsub__(a, b): return Fr.c(b) - a
    def __mul__(a, b):
        b = Fr.c(b); return Fr(a.n * b.n, a.d * b.d)
    __rmul__ = __mul__
    def __truediv__(a, b):
        b = Fr.c(b); assert b.n != 0; return Fr(a.n * b.d, a.d * b.n)
    def __rtruediv__(a, b): return Fr.c(b) / a
    def __pow__(a, k):
        if k >= 0: return Fr(a.n ** k, a.d ** k)
        return Fr(a.d ** (-k), a.n ** (-k))
    def __eq__(a, b):
        b = Fr.c(b); return (a.n * b.d - b.n * a.d) == 0
    def __ne__(a, b): return not a.__eq__(b)
    def as_expr(a): return sp.sympify(a.n.as_expr()) / sp.sympify(a.d.as_expr())

K = lambda x: Fr(R(x))
s, t = Fr(sR), Fr(tR)

# e-monomials: represented as dict {sorted tuple of positive indices: coeff in K}; e_0 = 1
def e(n):
    if n < 0: return {}
    if n == 0: return {(): K(1)}
    return {(n,): K(1)}

def emul(A, B):
    out = defaultdict(lambda: K(0))
    for a, x in A.items():
        for b, y in B.items():
            key = tuple(sorted(a + b)); out[key] = out[key] + x * y
    return dict(out)

def escale(c, A):
    return {k: c * v for k, v in A.items()}

def eadd(*As):
    out = defaultdict(lambda: K(0))
    for A in As:
        for k, v in A.items(): out[k] = out[k] + v
    return {k: v for k, v in out.items() if v != 0}

def poch(a, n):
    r = 1
    for i in range(n): r *= (1 - a * t**i)
    return r

def alpha(j):
    if j < 0: return 0
    r = 1
    for i in range(1, j + 1): r *= (s - t**i) / (1 - t**i)
    return r

def make_c(lesion=False):
    def c(n, j):
        if j < 0 or j > n: return 0
        sh = 1 if lesion else 0
        return poch(s, n - j) / poch(t, n - j) * (alpha(j) - s * t**(n - j + sh) * alpha(j - 1))
    return c

def make_F(c):
    def F(n, w):
        return sum(c(n, j) * w**j for j in range(n + 1))
    return F

def iszero(x):
    if isinstance(x, dict):
        return all(v == 0 for v in x.values())
    return Fr.c(x) == 0

ok = True
def report(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name); sys.stdout.flush()
    ok = ok and cond

c = make_c(); F = make_F(c)
cL = make_c(True); FL = make_F(cL)

# sanity: F_0 = 1, F_1 = (1-s)(1-tw)/(1-t)
for wv in [t**3, s/t, (s+2)/(t-3)]:
    report('F_0 = 1', iszero(F(0, wv) - 1))
    report('F_1 = (1-s)(1-tw)/(1-t)', iszero(F(1, wv) - (1 - s) * (1 - t * wv) / (1 - t)))

# Arm A
for n in range(1, 13):
    report(f'A: F_{n}(t^-{n}) = 0', iszero(F(n, t**(-n))))
    report(f'A control: F_{n}(t^-{n+1}) != 0', not iszero(F(n, t**(-(n + 1)))))

# Arm B
def coeff_formula(k, p, Ff):
    return eadd(*[escale(s**b * Ff(k - b, t**(p - b)), emul(e(b), e(p + k - b))) for b in range(0, k + 1)])

def coeff_gf(k, p, cc):
    # [z^p] of sum_{j+b<=k} s^b t^{-kj} c(k-b,j) e_b z^{b-k} E(t^j z); E(y) = sum e_n y^n
    tot = {}
    for b in range(0, k + 1):
        for j in range(0, k - b + 1):
            n = p + k - b  # need z^{b-k} * z^n = z^p
            if n < 0: continue
            tot = eadd(tot, escale(s**b * t**(-k * j) * cc(k - b, j) * t**(j * n), emul(e(b), e(n))))
    return tot

for k in range(1, 9):
    for N in range(1, k + 1):
        report(f'B formula: [z^-{N}]T_{k} = 0', iszero(coeff_formula(k, -N, F)))
        report(f'B gf:      [z^-{N}]T_{k} = 0', iszero(coeff_gf(k, -N, c)))
    report(f'B +control: [z^0]T_{k} = e_{k}', iszero(eadd(coeff_gf(k, 0, c), escale(K(-1), e(k)))))
    report(f'B -control: [z^1]T_{k} != 0', not iszero(coeff_gf(k, 1, c)))

# Arm C: lesion must fire
fired_i = sum(1 for n in range(1, 7) if not iszero(FL(n, t**(-n))))
fired_ii = sum(1 for k in range(1, 6) for N in range(1, k + 1) if not iszero(coeff_formula(k, -N, FL)))
tot_ii = sum(k for k in range(1, 6))
print(f'C lesion fires: (i) {fired_i}/6, (ii) {fired_ii}/{tot_ii}')
report('C lesion (i) fires on every n', fired_i == 6)
report('C lesion (ii) fires somewhere', fired_ii > 0)

# Arm D: Pieri rule vs subset operator, exact rationals
sv, tv = sp.Rational(3, 7), sp.Rational(5, 11)
def Ek_apply(k, r, m):
    X = sp.symbols(f'x1:{m+1}')
    def er(vars_):
        return sum(sp.Mul(*cmb) for cmb in itertools.combinations(vars_, r)) if r > 0 else sp.Integer(1)
    tot = 0
    for A in itertools.combinations(range(m), k):
        cA = 1
        for i in A:
            for jj in range(m):
                if jj not in A:
                    cA *= (X[i] - tv * X[jj]) / (X[i] - X[jj])
        xA = sp.Mul(*[X[i] for i in A])
        newv = [sv * X[i] if i in A else X[i] for i in range(m)]
        tot += cA * xA * er(newv)
    return X, sp.cancel(sp.together(tot))

def esym(X, n):
    if n < 0: return 0
    if n == 0: return 1
    return sum(sp.Mul(*cmb) for cmb in itertools.combinations(X, n))

nD = 0; okD = 0
for m in range(1, 5):
    for k in range(1, min(3, m) + 1):
        for r in range(0, 4):
            X, lhs = Ek_apply(k, r, m)
            def ev(x):
                return sp.sympify(x.as_expr()).subs({sS: sv, tS: tv})
            rhs = sum(ev(s**b * F(k - b, t**(r - b))) * esym(X, b) * esym(X, r + k - b)
                      for b in range(0, k + 1))
            nD += 1
            if sp.expand(lhs - rhs) == 0: okD += 1
print(f'D Pieri vs subset operator: {okD}/{nD}')
report('D transcription sanity', okD == nD)

print('ALL OK' if ok else 'SOME FAILED')

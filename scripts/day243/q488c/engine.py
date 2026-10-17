"""Rick's independent engine for Clio Q484 (Motegi arXiv:2610.15344 type B/C weight functions).

Route (deliberately different from Clio's): F_{N,n,r}(u|w|I) is computed from the printed
symmetrisation formula (Motherfunction, Thm explicitFformula) as an honest Laurent polynomial in
SYMBOLIC u (Vandermonde-cleared antisymmetrisation, exact division), and only then specialised at
u = w_M padded with 1's (n = N).  No use of the stability relation.
Normalisations E^C (typeCnormalizationfactor) and E^B (typeBnormalizationfactor) as printed.
"""
import itertools, sympy as sp
from functools import lru_cache

q = sp.Symbol('q')

def subsets(N):
    out = []
    for s in range(N + 1):
        out += [tuple(c) for c in itertools.combinations(range(1, N + 1), s)]
    return out

def geq(M, I):  # Motegi's order
    s, r = len(M), len(I)
    return s >= r and all(M[s - r + j] >= I[j] for j in range(r))

def summand(u, w, I, N, n, beta, gamma, qq, plant=None):
    r = len(I)
    t = sp.Integer(1)
    for j in range(r):
        ij = I[j]
        f = (1 - qq) * (u[j]**2 - 1)
        if plant == 'P2':
            f = (1 - qq) * u[j]**2
        for k in range(ij + 1, N + 1): f *= (u[j] - qq * w[k - 1])
        for k in range(1, ij): f *= (u[j] - w[k - 1])
        t *= f
    for j in range(r):
        for k in range(j + 1, r):
            t *= (qq * u[j] - u[k]) * (u[j] * u[k] - 1)   # /(u_j-u_k) handled by Vandermonde
    for j in range(r, n):
        for k in range(1, N + 1): t *= (u[j] - qq * w[k - 1])
    for j in range(r, n):
        for k in range(r):
            t *= (qq * u[j] - u[k]) * (u[j] * u[k] - 1)   # /(u_j-u_k), j>k
    for j in range(r, n):
        for k in range(j + 1, n):
            t *= (u[j] * u[k] - (1 if plant == 'P1' else qq))
    for j in range(r, n):
        t *= (beta * u[j] - gamma)
    return t

def denom_pairs(n, r):
    # pairs (a,b) meaning a factor 1/(u_a-u_b) in the unpermuted summand
    P = [(j, k) for j in range(r) for k in range(j + 1, r)]
    P += [(j, k) for j in range(r, n) for k in range(r)]
    return P

@lru_cache(maxsize=None)
def F_poly(N, n, I, beta, gamma, qq, plant=None):
    """returns (polynomial P(u), so that F = prod u^{-(N+n)} * P / (n-r)!) as sympy expr in u"""
    u = sp.symbols('u1:%d' % (n + 1))
    w = sp.symbols('w1:%d' % (N + 1))
    r = len(I)
    V = sp.Integer(1)
    for a in range(n):
        for b in range(a + 1, n): V *= (u[a] - u[b])
    tot = sp.Integer(0)
    for perm in itertools.permutations(range(n)):
        us = [u[p] for p in perm]
        t = summand(us, w, I, N, n, beta, gamma, qq, plant)
        # multiply by V / prod_{pairs}(us_a-us_b): the complementary pairs times sign
        pairs = set(denom_pairs(n, r))
        comp = sp.Integer(1)
        for a in range(n):
            for b in range(a + 1, n):
                if (a, b) in pairs: continue
                if (b, a) in pairs: comp *= -1; continue
                comp *= (us[a] - us[b])
        # V(u) vs V(us): V(us) = sign(perm) V(u); t*comp/V(us)... we accumulate t*comp*sign
        sgn = sp.combinatorics.Permutation(list(perm)).signature()
        tot += sgn * t * comp
    P, rem = sp.div(sp.Poly(sp.expand(tot), *u), sp.Poly(V, *u))
    assert rem.is_zero
    return P.as_expr() / sp.factorial(n - r), u, w

def W(X, N, I, M, qq=q, plant=None, normalise=True):
    """\\widetilde{WX}_{N,N,|I|}(w_M|w|I) (normalise=True) or WX (False). n=N."""
    n = N
    beta, gamma = (1, 0) if X == 'C' else (1, -1)
    P, u, w = F_poly(N, n, I, beta, gamma, qq, plant)
    uv = [w[m - 1] for m in M] + [sp.Integer(1)] * (n - len(M))
    sub = dict(zip(u, uv))
    val = P.subs(sub) * sp.prod([x**(-(N + n)) for x in uv])
    if X == 'B':
        val *= sp.prod([1 / (1 + 1 / x) for x in uv])
    if normalise:
        r = len(I)
        E = (1 - qq)**(-(n - r)) * sp.prod([1 - qq * a / b for a in uv for b in uv]) \
            * sp.prod([1 - qq * a * b for a in uv for b in w])
        E /= sp.prod([1 - qq * uv[j] * uv[k] for j in range(n) for k in range(j + 1, n)])
        if X == 'C':
            E /= sp.prod([1 - qq * w[j] * w[k] for j in range(N) for k in range(j, N)])
        else:
            E *= sp.prod([1 - qq * a for a in uv])
            E /= sp.prod([1 - qq * a**2 for a in uv])
            E /= sp.prod([1 - qq * w[j] * w[k] for j in range(N) for k in range(j + 1, N)])
            E /= sp.prod([1 - qq * b for b in w])
        val = val / E
    assert not val.atoms(sp.Float), "float leaked"
    return sp.factor(sp.cancel(val)), w

def structure_constants(X, N, qq=q, plant=None, normalise=True, lam_list=None):
    S = subsets(N)
    A = {}
    for M in S:
        for I in S:
            A[M, I], w = W(X, N, I, M, qq, plant, normalise)
    # solve A c = v for each (lam,mu) by forward substitution in a linear extension
    order = sorted(S, key=lambda T: (len(T), sum(T)))  # linear extension of >= ? verify
    for a in S:
        for b in S:
            if geq(a, b) and a != b: assert order.index(b) < order.index(a)
    C = {}
    pairs = lam_list if lam_list else [(l, m) for l in S for m in S]
    for (l, m) in pairs:
        c = {}
        for nu in order:
            rhs = A[nu, l] * A[nu, m] - sum(c[k] * A[nu, k] for k in c)
            c[nu] = sp.factor(sp.cancel(rhs / A[nu, nu]))
        C[l, m] = c
    return A, C, w

def is_laurent(expr, w):
    expr = sp.cancel(sp.together(expr))
    if expr == 0: return True
    num, den = sp.fraction(expr)
    den = sp.factor(den)
    # den must be (const in q) * monomial in w
    for f, e in sp.factor_list(den)[1]:
        if f.free_symbols & set(w):
            if not (sp.Poly(f, *w).is_monomial): return False
    return True

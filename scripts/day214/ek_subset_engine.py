"""Day 214: independent engine for e_k * F via the subset formula (207b A_k + K_k):
   E_k F = sum_{|A|=k} prod_{i in A, j notin A} (x_i - t x_j)/(x_i - x_j) * X_A * F(X_{A^c}, s X_A).
Exact over QQ (s,t numeric rationals or symbolic). Used to check DS (S),(L),(V) for all lengths."""
import sympy as sp, itertools, sys
def setup(m):
    xs = sp.symbols(f'x1:{m+1}'); return xs
def Ek(F, k, xs, s, t):
    m = len(xs); V = sp.prod([xs[i]-xs[j] for i in range(m) for j in range(i+1, m)])
    tot = 0
    for A in itertools.combinations(range(m), k):
        Ac = [j for j in range(m) if j not in A]
        sign = (-1)**sum(1 for i in A for j in Ac if i > j)
        w = sp.prod([xs[i]-xs[j] for i in A for j in A if i < j]) * sp.prod([xs[i]-xs[j] for i in Ac for j in Ac if i < j])
        num = sp.prod([xs[i]-t*xs[j] for i in A for j in Ac]) * sp.prod([xs[i] for i in A])
        sub = {xs[i]: s*xs[i] for i in A}
        tot += sign*w*num*sp.sympify(F).subs(sub, simultaneous=True)
    q, r = sp.div(sp.Poly(sp.expand(tot), *xs), sp.Poly(V, *xs))
    assert r.is_zero
    return q.as_expr()
def e(r, xs): return sum(sp.prod(c) for c in itertools.combinations(xs, r)) if r <= len(xs) else 0
def transpose(k):
    return tuple(sum(1 for p in k if p > i) for i in range(k[0])) if k else ()
def e_expand(P, xs):
    P = sp.Poly(sp.expand(P), *xs); out = {}
    while not P.is_zero:
        mon, c = max(P.terms(), key=lambda mc: mc[0])  # lex-leading monomial; exponents weakly decreasing for symmetric P
        kap = tuple(a for a in mon if a > 0); nu = transpose(kap)
        out[nu] = out.get(nu, 0) + c
        P = P - sp.Poly(c*sp.prod([e(p, xs) for p in nu]), *xs)
    return out
def dominates(a, b):
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0; sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True
def nstat(l): return sum(i*p for i, p in enumerate(l))
def parts(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield (); return
    for p in range(min(n, mx), 0, -1):
        for r in parts(n-p, p): yield (p,)+r
def e_lam_star(lam, xs, s, t, order=None):
    F = sp.Integer(1)
    for k in (reversed(lam) if order is None else order): F = Ek(F, k, xs, s, t)
    return F

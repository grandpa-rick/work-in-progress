"""Check Macdonald III (3.2): e_k P_rho(x;t) = sum_{kappa/rho vert k-strip} prod_v [m_v(kappa); r_v]_t P_kappa(x;t),
P via Macdonald III (2.2) symmetrization over S_m/S_m^lambda, in m = |kappa| variables."""
import sys, itertools, sympy as sp
from collections import Counter
sys.path.insert(0, '/home/agent/projects/scripts/day214'); sys.path.insert(0, '.')
from ek_subset_engine import e, parts
def qbin(N, r, t):
    if r < 0 or r > N: return 0
    return sp.cancel(sp.prod([1 - t**(N - i) for i in range(r)]) / sp.prod([1 - t**(i + 1) for i in range(r)]))
def P(lam, xs, t):
    m = len(xs); lam = list(lam) + [0] * (m - len(lam))
    V = sp.prod([xs[i] - xs[j] for i in range(m) for j in range(i + 1, m)]); tot = 0
    for w in set(itertools.permutations(lam)):  # w = exponent vector = w(lam); cosets <-> distinct rearrangements
        # term: x^w prod_{w_i>w_j} (x_i - t x_j)/(x_i - x_j); multiply by V
        num = sp.prod([xs[i] ** w[i] for i in range(m)])
        for i in range(m):
            for j in range(m):
                if w[i] > w[j]: num *= (xs[i] - t * xs[j]) * (1 if i < j else -1)
                elif w[i] == w[j] and i < j: num *= (xs[i] - xs[j])
                elif w[i] < w[j] and i < j: pass
        tot += num
    q, r = sp.div(sp.Poly(sp.expand(tot), *xs), sp.Poly(V, *xs)); assert r.is_zero
    return q
def expand_P(G, xs, t, n):
    out = {}; G = sp.Poly(G, *xs)
    while not G.is_zero:
        mon, c = max(G.terms(), key=lambda mc: mc[0]); kap = tuple(a for a in mon if a > 0)
        out[kap] = sp.factor(c); G = G - c * P(kap, xs, t)
    return out
t = sp.Rational(sys.argv[2]) if len(sys.argv) > 2 else sp.Symbol('t'); N = int(sys.argv[1]); bad = 0
for n in range(1, N + 1):
    xs = sp.symbols(f'x1:{n+1}')
    for k in range(1, n + 1):
        for rho in parts(n - k):
            G = sp.expand(e(k, xs) * P(rho, xs, t).as_expr()); got = expand_P(G, xs, t, n)
            r0 = list(rho) + [0] * k; pr = {}
            for A in itertools.combinations(range(len(r0)), k):
                kap = r0[:]
                for a in A: kap[a] += 1
                if any(kap[i] < kap[i + 1] for i in range(len(kap) - 1)): continue
                kp = tuple(x for x in kap if x > 0); mv = Counter(kp); rv = Counter(r0[a] + 1 for a in A)
                pr[kp] = sp.prod([qbin(mv[v], rv[v], t) for v in rv])
            ok = set(got) == set(pr) and all(sp.simplify(got[x] - pr[x]) == 0 for x in got)
            if not ok: bad += 1; print('FAIL', rho, k, got, pr)
    print('n', n, 'bad', bad, flush=True)

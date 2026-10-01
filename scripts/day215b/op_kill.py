"""Day 215 PROVE kill test: L_k = E_k mod s on b_mu = s^{n(mu)} e_mu.
Prediction (Theorem H, combinatorial form): with rho=mu', kappa=nu',
 L_k b_mu = sum_{kappa/rho vertical k-strip} t^{#{i in A,j notin A: kappa_i<kappa_j}} prod_v [m_v(kappa); r_v]_t b_nu
 where r_v = # parts of kappa equal to v that are decremented."""
import sys, itertools, sympy as sp
sys.path.insert(0, '/home/agent/projects/scripts/day214')
from ek_subset_engine import Ek, e, e_expand, transpose, nstat, parts
s, t = sp.symbols('s t')
def qbin(N, r, t):
    if r < 0 or r > N: return 0
    num = sp.prod([1 - t**(N - i) for i in range(r)]); den = sp.prod([1 - t**(i + 1) for i in range(r)])
    return sp.cancel(num / den)
def predicted(mu, k, tt):
    rho = transpose(mu) if mu else (); out = {}
    # kappa = rho + 1 in k distinct rows (vertical strip): choose counts per value level of kappa
    m = len(rho) + k
    r0 = list(rho) + [0] * k
    for A in itertools.combinations(range(m), k):
        kap = r0[:]
        for a in A: kap[a] += 1
        if any(kap[i] < kap[i + 1] for i in range(m - 1)): continue
        kap = tuple(x for x in kap if x > 0)
        # coefficient computed from kappa and r_v
        from collections import Counter
        mv = Counter(kap); rv = Counter(r0[a] + 1 for a in A)
        lt = sum(rv[v] * (mv[w] - rv[w]) for v in rv for w in mv if v < w)
        c = tt**lt * sp.prod([qbin(mv[v], rv[v], tt) for v in rv])
        nu = transpose(kap); out[nu] = sp.expand(out.get(nu, 0) + c)
    return out
N = int(sys.argv[1]); tt = sp.Rational(sys.argv[2]) if len(sys.argv) > 2 else t
bad = 0; tot = 0
for n in range(0, N):
    for mu in parts(n):
        for k in range(1, N - n + 1):
            m = n + k; xs = sp.symbols(f'x1:{m+1}')
            F = s**nstat(mu) * sp.prod([e(p, xs) for p in mu])
            G = Ek(F, k, xs, s, tt)
            ex = e_expand(G, xs); got = {}
            for nu, c in ex.items():
                c = sp.factor(c); cc = sp.cancel(c / s**nstat(nu))
                num, den = sp.fraction(cc)
                assert den.subs(s, 0) != 0, ('neg s power', mu, k, nu, cc)
                v = sp.simplify(cc.subs(s, 0))
                if v != 0: got[nu] = sp.expand(v)
            pr = predicted(mu, k, tt); tot += 1
            ok = set(got) == set(pr) and all(sp.simplify(got[x] - pr[x]) == 0 for x in got)
            if not ok: bad += 1; print('FAIL', mu, k, got, pr)
    print('done n+k<=', N, 'n=', n, 'bad', bad, 'tot', tot, flush=True)

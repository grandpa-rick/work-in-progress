"""Day 215: tabulate d_{lam mu}(t) = [s^{n(mu)}] c_{lam mu}(s,t), where e_lam^* = sum_mu c_{lam mu} e_mu.

Method (exact Fractions throughout):
 (1) one-step matrices: for each k and partition mu of n-k, a_{k mu nu}(s,t) = [e_nu] E_k e_mu computed by
     point evaluation of the subset formula (x-points exact rationals, n variables, as in day214/ds_pointeval.py),
     then exact interpolation in s (deg <= n-k) and in t (deg <= k(n-k)), each with one extra verification point.
 (2) rescaled basis b_mu = s^{n(mu)} e_mu: N_k[mu->nu](t) = [s^{n(nu)-n(mu)}] a_{k mu nu}; we CHECK that
     no lower power of s occurs (integrality), so reduction mod s is multiplicative and
     d_lam = N_{lam_1} N_{lam_2} ... N_{lam_l} (1).
 (3) independent cross-check (n<=5, a few t values): direct interpolation in s of c_{lam mu}(s,t0) from
     day214 coeffs(), lowest coefficient compared to d_{lam mu}(t0).
"""
import sys, itertools, math, random, pickle
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/scripts/day214')
from ds_pointeval import parts, dominates, nstat, esym, solve, coeffs

def poly_interp(xs, ys):
    """exact coefficient list of the interpolating polynomial"""
    return solve([[x**j for j in range(len(xs))] for x in xs], ys)

def peval(c, x): return sum(ci * x**i for i, ci in enumerate(c))

def trim(c):
    c = list(c)
    while c and c[-1] == 0: c.pop()
    return c

class Basis:
    """for fixed number of variables n: inverse of e_nu point-evaluation matrix, reused"""
    def __init__(self, n, rng):
        self.n = n; self.P = list(parts(n))
        self.pts = [tuple(Fr(rng.randint(-10**4, 10**4), rng.randint(1, 10**3)) for _ in range(n)) for _ in self.P]
        self.M = [[math.prod([esym(p, x) for p in nu]) for nu in self.P] for x in self.pts]
        self.chk = tuple(Fr(rng.randint(-10**4, 10**4), rng.randint(1, 10**3)) for _ in range(n))
    def expand(self, f):
        c = solve(self.M, [f(x) for x in self.pts])
        x = self.chk
        assert sum(ci*math.prod([esym(p, x) for p in nu]) for ci, nu in zip(c, self.P)) == f(x)
        return dict(zip(self.P, c))

def op(k, mu, x, s, t):
    m = len(x); tot = Fr(0)
    for A in itertools.combinations(range(m), k):
        As = set(A); w = Fr(1)
        for i in A:
            w *= x[i]
            for j in range(m):
                if j not in As: w *= (x[i]-t*x[j])/(x[i]-x[j])
        y = [s*x[i] if i in As else x[i] for i in range(m)]
        tot += w*math.prod([esym(p, y) for p in mu])
    return tot

def one_step(B, k, mu):
    """returns {nu: N_k[mu->nu](t) as coeff list}, after checking integrality in b-basis"""
    n = B.n; ds, dt = n-k, k*(n-k)
    sv = [Fr(i+2, 3) for i in range(ds+2)]
    tv = [Fr(-(i+1), 5) if i % 2 else Fr(i+1, 4) for i in range(dt+2)]
    data = {(a, b): B.expand(lambda x: op(k, mu, x, a, b)) for a in sv for b in tv}
    out = {}
    for nu in B.P:
        # interpolate in t for each s, then in s for each t-coefficient
        tpolys = []
        for a in sv:
            c = poly_interp(tv[:-1], [data[(a, b)][nu] for b in tv[:-1]])
            assert peval(c, tv[-1]) == data[(a, tv[-1])][nu], 't-degree bound fail'
            tpolys.append(c)
        # coefficient matrix C[i][j] = [s^i t^j]
        C = [[None]*(dt+1) for _ in range(ds+1)]
        for j in range(dt+1):
            c = poly_interp(sv[:-1], [tp[j] for tp in tpolys[:-1]])
            assert peval(c, sv[-1]) == tpolys[-1][j], 's-degree bound fail'
            for i in range(ds+1): C[i][j] = c[i]
        off = nstat(nu) - nstat(mu)
        for i in range(ds+1):
            if i < off and any(C[i]): raise AssertionError(f'integrality FAIL k={k} mu={mu} nu={nu} s^{i}')
        if 0 <= off <= ds and any(C[off]): out[nu] = trim(C[off])
    return out

def pmul(a, b):
    if not a or not b: return []
    r = [Fr(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
def padd(a, b):
    r = [Fr(0)]*max(len(a), len(b))
    for i, x in enumerate(a): r[i] += x
    for i, x in enumerate(b): r[i] += x
    return trim(r)

if __name__ == '__main__':
    n = int(sys.argv[1]); rng = random.Random(215 + n)
    # need one-step matrices with n' = |mu|+k for all n' <= n ... but ALL steps use n variables (as day214 star()).
    B = Basis(n, rng)
    # intermediate states live in partitions of m<=n, but evaluated in n variables; one_step for size m uses n vars
    # generalise: one_step with nu ranging over partitions of k+|mu| -> build per-size bases in n variables.
    bases = {}
    def basis_for(m):
        if m not in bases:
            b = Basis.__new__(Basis); b.n = n; b.P = list(parts(m))
            b.pts = [tuple(Fr(rng.randint(-10**4, 10**4), rng.randint(1, 10**3)) for _ in range(n)) for _ in b.P]
            b.M = [[math.prod([esym(p, x) for p in nu]) for nu in b.P] for x in b.pts]
            b.chk = tuple(Fr(rng.randint(-10**4, 10**4), rng.randint(1, 10**3)) for _ in range(n))
            bases[m] = b
        return bases[m]
    def one_step_m(k, mu):
        m = k + sum(mu); b = basis_for(m)
        # degree bounds: s-degree <= |mu|; t-degree <= k(n-k) (n variables)
        ds, dt = sum(mu), k*(n-k)
        sv = [Fr(i+2, 3) for i in range(ds+2)]
        tv = [Fr(-(i+1), 5) if i % 2 else Fr(i+1, 4) for i in range(dt+2)]
        data = {(a, c): b.expand(lambda x: op(k, mu, x, a, c)) for a in sv for c in tv}
        out = {}
        for nu in b.P:
            tpolys = []
            for a in sv:
                cc = poly_interp(tv[:-1], [data[(a, c)][nu] for c in tv[:-1]])
                assert peval(cc, tv[-1]) == data[(a, tv[-1])][nu]
                tpolys.append(cc)
            C = [[None]*(dt+1) for _ in range(ds+1)]
            for j in range(dt+1):
                cc = poly_interp(sv[:-1], [tp[j] for tp in tpolys[:-1]])
                assert peval(cc, sv[-1]) == tpolys[-1][j]
                for i in range(ds+1): C[i][j] = cc[i]
            off = nstat(nu) - nstat(mu)
            for i in range(min(off, ds+1)):
                assert not any(C[i]), f'integrality FAIL k={k} mu={mu} nu={nu} s^{i}'
            if 0 <= off <= ds and any(C[off]): out[nu] = trim(C[off])
        return out
    cache = {}
    def N(k, mu):
        if (k, mu) not in cache: cache[(k, mu)] = one_step_m(k, mu)
        return cache[(k, mu)]
    D = {}
    for lam in parts(n):
        state = {(): [Fr(1)]}
        for k in reversed(lam):
            new = {}
            for mu, p in state.items():
                for nu, q in N(k, mu).items():
                    new[nu] = padd(new.get(nu, []), pmul(p, q))
            state = {nu: p for nu, p in new.items() if p}
        D[lam] = state
        print(lam, {mu: [str(x) for x in p] for mu, p in state.items()}, flush=True)
    pickle.dump(D, open(f'/home/agent/projects/scripts/day215/d_n{n}.pkl', 'wb'))
    pickle.dump(cache, open(f'/home/agent/projects/scripts/day215/N_n{n}.pkl', 'wb'))
    # cross-check vs direct computation (day214 coeffs, n variables) at t0 values
    if len(sys.argv) > 2:
        for t0 in [Fr(2), Fr(-2, 7)]:
            for lam in parts(n):
                Ds = nstat(lam) + n + 1
                sv = [Fr(i+2, 3) for i in range(Ds+2)]
                data = [coeffs(lam, x, t0, rng) for x in sv]
                for mu in parts(n):
                    c = poly_interp(sv[:-1], [d.get(mu, Fr(0)) for d in data[:-1]])
                    assert peval(c, sv[-1]) == data[-1].get(mu, Fr(0))
                    v = next((i for i, a in enumerate(c) if a != 0), None)
                    if dominates(mu, lam):
                        assert v >= nstat(mu) and c[nstat(mu)] == peval(D[lam][mu], t0), (lam, mu, t0)
                    else:
                        assert v is None
            print('direct cross-check OK at t0 =', t0, flush=True)

"""Day 214: fast exact check of DS via point evaluation of the subset formula (0.1).
E_k F (x) = sum_{|A|=k} prod_{i in A, j notin A} (x_i - t x_j)/(x_i - x_j) * X_A * F(x with x_a -> s x_a, a in A).
e_lam^{(q,t)}(x) evaluated recursively (exact Fractions); e-coefficients by solving a linear system over Q.
Checks (S),(L),(V) and full-up-set support at a fixed exact (s,t). Nonzero value at a rational point => nonzero coefficient."""
from fractions import Fraction as Fr
import itertools, random, sys
from functools import lru_cache
def parts(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield (); return
    for p in range(min(n, mx), 0, -1):
        for r in parts(n-p, p): yield (p,)+r
def dominates(a, b):
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0; sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True
def nstat(l): return sum(i*p for i, p in enumerate(l))
def esym(r, x):
    E = [Fr(1)] + [Fr(0)]*len(x)
    for v in x:
        for j in range(len(x), 0, -1): E[j] += v*E[j-1]
    return E[r] if r <= len(x) else Fr(0)
def star(ks, x, s, t):
    # E_{ks[0]} E_{ks[1]} ... (1) at point x
    if not ks: return Fr(1)
    k, m = ks[0], len(x); tot = Fr(0)
    for A in itertools.combinations(range(m), k):
        As = set(A); w = Fr(1)
        for i in A:
            w *= x[i]
            for j in range(m):
                if j not in As: w *= (x[i]-t*x[j])/(x[i]-x[j])
        y = tuple(s*x[i] if i in As else x[i] for i in range(m))
        tot += w*star(ks[1:], y, s, t)
    return tot
def solve(M, b):
    n = len(M); A = [row[:] + [bb] for row, bb in zip(M, b)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]/A[c][c]; A[r] = [a - f*bq for a, bq in zip(A[r], A[c])]
    return [A[i][n]/A[i][i] for i in range(n)]
def coeffs(lam, s, t, rng):
    n = sum(lam); P = list(parts(n)); N = len(P)
    pts = [tuple(Fr(rng.randint(-10**6, 10**6), rng.randint(1, 10**4)) for _ in range(n)) for _ in range(N)]
    M = [[__import__('math').prod([esym(p, x) for p in nu]) for nu in P] for x in pts]
    b = [star(list(lam), x, s, t) for x in pts]
    c = solve(M, b)
    # extra-point consistency check
    x = tuple(Fr(rng.randint(-10**6, 10**6), rng.randint(1, 10**4)) for _ in range(n))
    assert sum(ci*__import__('math').prod([esym(p, x) for p in nu]) for ci, nu in zip(c, P)) == star(list(lam), x, s, t)
    return {nu: ci for nu, ci in zip(P, c) if ci != 0}
if __name__ == '__main__':
    n0, n1 = int(sys.argv[1]), int(sys.argv[2]); rng = random.Random(214)
    S, T = Fr(3, 7), Fr(-5, 11); ok = True; nfull = ntot = 0
    for n in range(n0, n1+1):
        for lam in parts(n):
            out = coeffs(lam, S, T, rng)
            supp = all(dominates(mu, lam) for mu in out)
            lead = out.get(lam, 0) == S**nstat(lam)
            V = coeffs(lam, Fr(1), T, rng) == {lam: 1}
            rev = coeffs(tuple(reversed(lam)), S, T, rng) == out  # operator order independence
            full = set(out) == {mu for mu in parts(n) if dominates(mu, lam)}
            ok &= supp and lead and V and rev; nfull += full; ntot += 1
            print(lam, 'supp', supp, 'lead', lead, 's=1', V, 'order-indep', rev, 'full-upset', full, flush=True)
    print(f'full-upset {nfull}/{ntot}'); print('ALL OK' if ok else 'FAIL')

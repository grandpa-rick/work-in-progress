# Wake 253 Q27: OG(m,2m+1) (type B_m, spinor O(1)=L_{omega_m}) Schubert varieties X_kappa (dimension-indexed, kappa strict, parts in [m]).
#  (1) H^OG_kappa(j) from minuscule SMT model  Omega(P_kappa, j+1) = #order-preserving P_kappa -> {0..j}  (P_kappa = shifted diagram)
#      CHECKED against Demazure character formula in type B_m (independent), and full OG against Weyl dimension formula.
#  (2) [O_{X_kappa}].O^1 : BR Cor 4.8 (KOG-tableaux, content {1..p}) implemented for general p; at p=1 compared with the TRUE K-class
#      computed by pulling back to G/B and solving with Demazure characters (as wake252 kclass.py, type B).
#  (3) zeros of H^OG_kappa vs prediction {-1,...,-r(kappa)}.
import sys, itertools, collections, random
from fractions import Fraction as Fr
from functools import lru_cache
from math import prod
import numpy as np
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/prove249d')
from orderpoly import shifted, boxless, longest_chain, ideals
from dem import word
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
RR = 2 if N <= 4 else 3
MODE = sys.argv[2] if len(sys.argv) > 2 else 'all'
def strict_parts(n): return [tuple(sorted(c, reverse=True)) for L in range(n + 1) for c in itertools.combinations(range(1, n + 1), L)]
def nm(t): return ''.join(map(str, t)) or '0'
# ---------- type B_n Demazure, doubled coordinates M = 2 mu ----------
def DiB(f, i, n):
    g = collections.defaultdict(int)
    for mu, a in f.items():
        if i < n: k = (mu[i - 1] - mu[i]) // 2; al = [0] * n; al[i - 1] = 2; al[i] = -2
        else: k = mu[n - 1]; al = [0] * n; al[n - 1] = 2       # alpha_n = e_n (short), <mu,alpha_n^vee> = 2 mu_n = M_n
        if k >= 0:
            for t in range(k + 1): g[tuple(mu[j] - t * al[j] for j in range(n))] += a
        elif k <= -2:
            for t in range(1, -k): g[tuple(mu[j] + t * al[j] for j in range(n))] -= a
    return {k: v for k, v in g.items() if v}
@lru_cache(None)
def wordk(k, n): return tuple(word(set(k), n))
def chi_dem(k, n, j):      # chi(X_kappa, O(j)) on OG, Demazure
    f = {tuple([j] * n): 1}
    for i in reversed(wordk(k, n)): f = DiB(f, i, n)
    return sum(f.values())
def weyl_B_spin(n, j):    # dim V(j omega_n), B_n ; rho = (n-1/2, ..., 1/2)
    lam = [Fr(j, 2)] * n; rho = [Fr(2 * (n - i) - 1, 2) for i in range(n)]
    num = den = Fr(1)
    for a in range(n):
        for b in range(a + 1, n):
            for s in (1, -1):    # e_a - e_b, e_a + e_b  (coroot = root for long)
                num *= (lam[a] + rho[a]) + s * -(lam[b] + rho[b]) if False else ((lam[a] + rho[a]) - s * (lam[b] + rho[b]))
                den *= (rho[a] - s * rho[b])
        num *= 2 * (lam[a] + rho[a]); den *= 2 * rho[a]   # short root e_a, coroot 2e_a
    return num / den
# ---------- order polynomial model ----------
@lru_cache(None)
def omega_vals(k, Jmax):     # Omega(P_k, t) for t=1..Jmax+1 i.e. H(j), j=0..Jmax
    P = shifted(list(k)); J = ideals(P, boxless); full = frozenset(P)
    sub = {I: [K for K in J if I <= K] for I in J}
    vec = {I: 1 for I in J}; out = [vec.get(full, 0)]
    for v in range(1, Jmax + 1):
        nv = collections.defaultdict(int)
        for I, c in vec.items():
            for K in sub[I]: nv[K] += c
        vec = nv; out.append(vec.get(full, 0))
    return tuple(out)
def H(k, x):
    v = omega_vals(k, sum(k)); n = len(v); s = Fr(0)
    for i in range(n):
        t = Fr(v[i])
        for l in range(n):
            if l != i: t *= Fr(x - l, i - l)
        s += t
    return s
def r(k): return longest_chain(shifted(list(k)), boxless) if k else 0
# ---------- BR KOG rule ----------
def boxes(l): return {(i, i + j) for i, li in enumerate(l) for j in range(li)}
def dual(k, n): return tuple(sorted(set(range(1, n + 1)) - set(k), reverse=True))
def is_rim(th): return not any(r2 > r1 and c2 > c1 for (r1, c1) in th for (r2, c2) in th)
def kog_count(th, p):
    if not is_rim(th): return 0
    th = sorted(th); n = 0
    for lab in itertools.product(range(1, p + 1), repeat=len(th)):
        if set(lab) != set(range(1, p + 1)): continue
        L = dict(zip(th, lab)); ok = True
        for (a, b), x in L.items():
            if (a, b + 1) in L and not x < L[(a, b + 1)]: ok = False; break
            if (a + 1, b) in L and not x < L[(a + 1, b)]: ok = False; break
            sw = [y for (c, d), y in L.items() if c > a and d < b]
            if sw and not (all(x <= y for y in sw) or all(x >= y for y in sw)): ok = False; break
        n += ok
    return n
def kog_pieri_dim(k, n, p=1):
    lam = dual(k, n); bl = boxes(lam); out = {}
    for nu in strict_parts(n):
        bn = boxes(nu)
        if bl <= bn and bl != bn:
            th = bn - bl; c = (-1) ** (len(th) - p) * kog_count(th, p)
            if c: out[dual(nu, n)] = c
    return out
# ---------- true K-class via G/B ----------
W0P = None
def chiGB(k, lam, n):
    f = {lam: 1}
    for i in [i for t in range(n - 1) for i in range(1, n - t)]: f = DiB(f, i, n)
    for i in reversed(wordk(k, n)): f = DiB(f, i, n)
    return sum(f.values())
def leq(a, b): return len(a) <= len(b) and all(a[i] <= b[i] for i in range(len(a)))
def true_kclass(k, n, LAMS, cache):
    def c(rr, l):
        if (rr, l) not in cache: cache[(rr, l)] = chiGB(rr, l, n)
        return cache[(rr, l)]
    rh = [q for q in strict_parts(n) if leq(q, k) and q != k]
    A = [[c(q, l) for q in rh] for l in LAMS]
    b = [c(k, l) - c(k, tuple(x - 1 for x in l)) for l in LAMS]
    An = np.array(A, dtype=float); assert np.linalg.matrix_rank(An) == len(rh), ('not separating', k)
    sol = [int(round(v)) for v in np.linalg.lstsq(An, np.array(b, dtype=float), rcond=None)[0]]
    assert all(sum(a * cc for a, cc in zip(row, sol)) == bb for row, bb in zip(A, b)), ('inexact', k)
    return {q: v for q, v in zip(rh, sol) if v}
def height(th, k):   # max |theta ∩ chain| over chains in P_kappa (longest chain in theta itself under boxless)
    return longest_chain(sorted(th), boxless) if th else 0
if __name__ == '__main__':
    n = N; KS = [k for k in strict_parts(n) if k]
    if MODE in ('all', 'model'):
        # full OG: Weyl dim vs Omega(staircase)
        st = tuple(range(n, 0, -1))
        okW = all(weyl_B_spin(n, j) == omega_vals(st, 8)[j] for j in range(9))
        print(f'm={n} full OG: Weyl dim V(j omega_m) == Omega(staircase,j+1) for j=0..8:', okW)
        okD = 0; jj = 3 if n <= 4 else 2
        for k in KS:
            d = [chi_dem(k, n, j) for j in range(jj + 1)]; o = list(omega_vals(k, max(sum(k), jj))[:jj + 1])
            okD += d == o
            if d != o: print('DEMAZURE MISMATCH', k, d, o)
        print(f'm={n}: Demazure(type B) == Omega model at j=0..{jj}: {okD}/{len(KS)}')
    if MODE in ('all', 'kclass'):
        random.seed(1); LAMS = []
        for _ in range(6 * 2 ** n):
            par = random.randint(0, 1); LAMS.append(tuple(2 * random.randint(-RR, RR) + par for _ in range(n)))
        cache = {}; ok = 0
        for k in KS:
            T = true_kclass(k, n, LAMS, cache); B = kog_pieri_dim(k, n, 1)
            ok += T == B
            if T != B: print('KCLASS MISMATCH', k, {nm(q): v for q, v in T.items()}, {nm(q): v for q, v in B.items()})
        print(f'm={n}: true K-class [O_X_kappa].O^1 (G/B Demazure solve) == BR KOG rule p=1: {ok}/{len(KS)}')
    if MODE in ('all', 'zeros'):
        hist = collections.Counter(); cons = 0; fails = []; rook = 0; hz = collections.Counter(); halfz = []
        for k in KS:
            K = kog_pieri_dim(k, n, 1)
            rook += all(not any(a == c or b == d for (a, b) in th for (c, d) in th if (a, b) != (c, d))
                        for th in (boxes(dual(q, n)) - boxes(dual(k, n)) for q in K))
            for q in K: hist[height(boxes(k) - boxes(q), k)] += 1
            pts = [Fr(t, 2) for t in range(-4 * n - 4, 2 * sum(k) + 4)]
            cons += all(H(k, x) - H(k, x - 1) == sum(c * H(q, x) for q, c in K.items()) for x in pts)
            R = sum(k) + 4
            zi = [j for j in range(-R, 0) if H(k, j) == 0]; zh = [Fr(2 * j + 1, 2) for j in range(-R, 0) if H(k, Fr(2 * j + 1, 2)) == 0]
            pred = list(range(-r(k), 0)); ok = zi == pred
            if not ok: fails.append((k, zi, pred))
            if zh: halfz.append((k, zh))
            print(nm(k), '|k|', sum(k), 'r', r(k), "r'", r(k[1:]), 'int zeros', zi, 'half-int zeros', [str(z) for z in zh],
                  'pred ok', ok, 'H(-r-1)', H(k, -r(k) - 1), 'supp', {nm(q): v for q, v in K.items()}, flush=True)
        print(f'm={n}: supp all rook strips {rook}/{len(KS)}; height (max |theta∩chain|) histogram over all (kappa,rho):', sorted(hist.items()))
        print(f'm={n}: b=H(j)-H(j-1) == sum_rho c_rho H_rho on half-int grid (KOG rule x Omega model): {cons}/{len(KS)}')
        print(f'm={n}: integer zero set == {{-1..-r(kappa)}} exactly: {len(KS)-len(fails)}/{len(KS)}; failures {fails}')
        print(f'm={n}: kappa with half-integer zeros: {len(halfz)}', halfz[:20])

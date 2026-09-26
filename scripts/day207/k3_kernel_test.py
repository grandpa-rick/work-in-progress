"""Day 207: |A|=3 parabolic kernel test for e_3(Y).e_r  (generalizes day206/parabolic_kernel_test.py).

Conventions (verbatim Day 198/205 pipeline, flint class day205/k3_fast_pipeline.AHA; s = q^{-1}):
  T_j F = t s_j F - (t-1) X_{j+1}(F - s_jF)/(X_j - X_{j+1}),  pi F = X_1 F(X_2..X_m, s X_1),
  Y_i = t^{m-i} T_{i-1}..T_1 pi T_{m-1}^{-1}..T_i^{-1}.
LHS  P  := e_k(Y).e_r = sum_{i1<..<ik} Y_i1..Y_ik e_r(X_1..X_m)   (k = 2 calibration, k = 3 test)
S       := sum_{w in W^J} T_w pi^k e_r, W^J = min coset reps of S_m/(S_k x S_{m-k}),
           T_w = T_{a1-1}..T_1 T_{a2-1}..T_2 ... T_{ak-1}..T_k  for 1<=a1<..<ak<=m.
K_conv  := sum_{|A|=k} F^{(A)} prod_{i in A, j notin A} kern(X_i,X_j), F = pi^k e_r
           = X_1..X_k e_r(X_{k+1}..X_m, sX_1..sX_k);  kern_T = (X_i - tX_j)/(X_i-X_j),
           kern_Tinv = (tX_i - X_j)/(X_i - X_j).
Pre-registered prediction (question file): sigma^(k) pi^k = t^{-C(k,2)} e_k(Y), i.e. P = t^{C(k,2)} S.
Test 2: (1-t)^3 R_(n,p,s) vs QJ(n,p,s) = [z^n w^p v^s] Q(z)Q(w)Q(v) prod_{i<j}(1-z_j/z_i)/(1-t z_j/z_i),
  R_alpha = sum_{a,b,c distinct} x_a^n x_b^p x_c^s prod_{i<j in (a,b,c)} kern_T(x_i,x_j) prod_{i in A, j notin A} kern_T.
"""
import sys, random, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA


def ekY_er(A_, k, r):
    """sum_{i1<..<ik} Y_i1..Y_ik e_r via nested suffix sums (Y's commute)."""
    m = A_.m
    G = {j: A_.e(r) for j in range(1, m + 2)}  # level 0: 'suffix sum' is just e_r
    for level in range(k):
        # H[i] = Y_i applied to (sum over tails starting after i) ; tail sum at level 0 = e_r
        Yi = {i: A_.Y(G[i + 1] if level > 0 else G[i], i) for i in range(1, m + 1)}
        # suffix sums: new G[i] = sum_{i' >= i} Yi[i']
        newG = {m + 1: 0 * A_.s}
        for i in range(m, 0, -1):
            newG[i] = newG[i + 1] + Yi[i]
        G = newG
    return G[1]


def coset_sym(A_, F, k):
    m = A_.m; tot = 0 * A_.s
    for a in itertools.combinations(range(1, m + 1), k):
        H = F
        for idx in range(k - 1, -1, -1):  # apply T_{ak-1}..T_k first, ..., T_{a1-1}..T_1 last
            for j in range(idx + 1, a[idx]):
                H = A_.T(H, j)
        tot = tot + H
    return tot


def evalp(A_, P, X, s, t):
    tot = Fr(0)
    for mon, c in P.to_dict().items():
        v = Fr(int(c))
        for i in range(A_.m): v *= X[i] ** int(mon[i])
        tot += v * s ** int(mon[A_.m]) * t ** int(mon[A_.m + 1])
    return tot


def e_vals(vals, r):
    if r < 0 or r > len(vals): return Fr(0)
    tot = Fr(0)
    for c in itertools.combinations(vals, r):
        p = Fr(1)
        for v in c: p *= v
        tot += p
    return tot


def kern(xi, xj, t, conv):
    return (xi - t * xj) / (xi - xj) if conv == 'T' else (t * xi - xj) / (xi - xj)


def kernel_sum(X, s, t, r, k, conv):
    m = len(X); tot = Fr(0)
    for A in itertools.combinations(range(m), k):
        w = Fr(1)
        for i in A:
            for j in range(m):
                if j not in A: w *= kern(X[i], X[j], t, conv)
        tail = [X[j] for j in range(m) if j not in A]
        F = e_vals(tail + [s * X[i] for i in A], r)
        for i in A: F *= X[i]
        tot += w * F
    return tot


def rand_pts(m, npts, seed):
    random.seed(seed)
    pts = []
    for _ in range(npts):
        X = [Fr(v, random.randint(1, 7)) for v in random.sample(range(2, 300), m)]
        while len(set(X)) < m:
            X = [Fr(v) for v in random.sample(range(2, 300), m)]
        s = Fr(random.randint(2, 9), random.randint(1, 5))
        t = Fr(random.randint(2, 9), random.randint(1, 5))
        if t == 1: t = Fr(3, 2)
        pts.append((X, s, t))
    return pts


def test1(m, r, k, npts=3):
    t0 = time.time()
    A_ = AHA(m)
    P = ekY_er(A_, k, r)
    F = A_.e(r)
    for _ in range(k): F = A_.pi(F)
    S = coset_sym(A_, F, k)
    res = {}
    for e in range(0, 5):  # exact polynomial identity P == t^e S ?
        res[f'P == t^{e} S (exact poly)'] = (P == A_.t ** e * S)
    pts = rand_pts(m, npts, 1000 * k + 10 * m + r)
    res['S == K_T (pts)'] = all(evalp(A_, S, X, s, t) == kernel_sum(X, s, t, r, k, 'T') for X, s, t in pts)
    for conv in ('T', 'Tinv'):
        rats = []
        for X, s, t in pts:
            R = kernel_sum(X, s, t, r, k, conv); L = evalp(A_, P, X, s, t)
            rats.append((L / R if R != 0 else None, t))
        res[f'LHS/K_{conv} (ratio, t)'] = rats
    res['P nonzero'] = not P.is_zero()
    res['time'] = round(time.time() - t0, 1)
    return res


# ---------------- Test 2: three-row functional vs 3-fold Jing product -------------
def qn_list(X, t, N):
    """q_n(X;t) = [z^n] prod (1 - t x z)/(1 - x z), n = 0..N."""
    c = [Fr(1)] + [Fr(0)] * N
    for x in X:
        # multiply by (1 - t x z) * sum (x z)^k
        g = [x ** n for n in range(N + 1)]
        g2 = [g[n] - (t * x * g[n - 1] if n >= 1 else 0) for n in range(N + 1)]
        c = [sum(c[i] * g2[n - i] for i in range(n + 1)) for n in range(N + 1)]
    return c


def QJ3(alpha, qn, t):
    n, p, s_ = alpha
    f = lambda kk: Fr(1) if kk == 0 else t ** kk - t ** (kk - 1)
    Q = lambda i: qn[i] if 0 <= i < len(qn) else Fr(0)
    tot = Fr(0); N = n + p + s_
    for k12 in range(N + 1):
        for k13 in range(N + 1 - k12):
            for k23 in range(N + 1):
                a, b, c = n + k12 + k13, p - k12 + k23, s_ - k13 - k23
                if b < 0 or c < 0: continue
                tot += f(k12) * f(k13) * f(k23) * Q(a) * Q(b) * Q(c)
    return tot


def R3(alpha, X, t):
    m = len(X); tot = Fr(0)
    for (a, b, c) in itertools.permutations(range(m), 3):
        A = (a, b, c)
        w = X[a] ** alpha[0] * X[b] ** alpha[1] * X[c] ** alpha[2]
        w *= kern(X[a], X[b], t, 'T') * kern(X[a], X[c], t, 'T') * kern(X[b], X[c], t, 'T')
        for i in A:
            for j in range(m):
                if j not in A: w *= kern(X[i], X[j], t, 'T')
        tot += w
    return tot


def test2(m, maxpart=4, minpart=0, npts=2):
    pts = rand_pts(m, npts, 77 + m)
    fails = []; n_ok = 0
    for alpha in itertools.product(range(minpart, maxpart + 1), repeat=3):
        ok = True
        for X, _, t in pts:
            qn = qn_list(X, t, sum(alpha))
            if (1 - t) ** 3 * R3(alpha, X, t) != QJ3(alpha, qn, t):
                ok = False; break
        if ok: n_ok += 1
        else: fails.append(alpha)
    return n_ok, fails


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('all', 't1'):
        cases = [(2, 5, 2), (2, 6, 2), (3, 6, 2), (3, 6, 0), (3, 6, 1), (3, 6, 3), (3, 7, 2), (3, 7, 3), (3, 7, 1), (3, 7, 4)]
        for k, m, r in cases:
            res = test1(m, r, k)
            print(f'k={k} m={m} r={r}:')
            for key, v in res.items():
                print(f'    {key}: {v}')
            sys.stdout.flush()
    if which in ('all', 't2'):
        for m in (4, 5, 6):
            t0 = time.time()
            n_ok, fails = test2(m, maxpart=4 if m < 6 else 3)
            print(f'Test2 m={m}: (1-t)^3 R_alpha == QJ3(alpha) for {n_ok} compositions (parts 0..{4 if m < 6 else 3}); fails={fails}  [{time.time()-t0:.1f}s]')
            sys.stdout.flush()

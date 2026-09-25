"""Day 206: test the |A|=2 parabolic HL kernel hypothesis for e_2(Y).e_r.

LHS  := e_2(Y).e_r = sum_{i<j} Y_i Y_j e_r(X_1..X_m), computed with the Day 205
        flint pipeline (day205/k3_fast_pipeline.AHA, conventions verbatim from
        day198/p2Y_er.build_action; s = q^{-1}).
RHS  := K^{(2)}(F), F = pi^2 e_r = X1 X2 e_r(X3..Xm, s X1, s X2),
        K^{(2)} F = sum_{|A|=2} F^{(A)} prod_{i in A, j notin A} kern(X_i, X_j)
        with kern = (X_i - t X_j)/(X_i - X_j)   [convention T ]
          or kern = (t X_i - X_j)/(X_i - X_j)   [convention T^{-1}]
Comparison: exact rational evaluation at random points (X, s, t); report ratio.
Also checks the operator identity  sum_{w in W^J} T_w F = K^{(2)}_T F  (flint, exact).
"""
import sys, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA

def e_vals(vals, r):
    if r < 0 or r > len(vals): return Fr(0)
    tot = Fr(0)
    for c in itertools.combinations(vals, r):
        p = Fr(1)
        for v in c: p *= v
        tot += p
    return tot

def F_pi2_er(X, s, r, heads):
    """(pi^2 e_r) with head variables (x_a, x_b) = values at positions heads, tail = rest."""
    a, b = heads
    tail = [X[k] for k in range(len(X)) if k not in heads]
    return X[a] * X[b] * e_vals(tail + [s * X[a], s * X[b]], r)

def kernel(X, s, t, r, conv, Ffun=F_pi2_er):
    m = len(X); tot = Fr(0)
    for A in itertools.combinations(range(m), 2):
        w = Fr(1)
        for i in A:
            for j in range(m):
                if j in A: continue
                if conv == 'T':
                    w *= (X[i] - t * X[j]) / (X[i] - X[j])
                else:
                    w *= (t * X[i] - X[j]) / (X[i] - X[j])
        tot += w * Ffun(X, s, r, A)
    return tot

def lhs_poly(A_, r):
    er = A_.e(r); m = A_.m
    Yv = {j: A_.Y(er, j) for j in range(1, m + 1)}
    tot = 0 * A_.s
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            tot = tot + A_.Y(Yv[j], i)
    return tot

def evalp(A_, P, X, s, t):
    tot = Fr(0)
    for mon, c in P.to_dict().items():
        v = Fr(int(c))
        for k in range(A_.m): v *= X[k] ** int(mon[k])
        v *= s ** int(mon[A_.m]) * t ** int(mon[A_.m + 1])
        tot += v
    return tot

def parabolic_sym_op(A_, F):
    """sum over minimal coset reps of S_m/(S_2 x S_{m-2}): T_{a-1}..T_1 T_{b-1}..T_2 F, 1<=a<b<=m."""
    m = A_.m; tot = 0 * A_.s
    for b in range(2, m + 1):
        G = F
        for k in range(2, b): G = A_.T(G, k)
        for a in range(1, b):
            H = G
            for k in range(1, a): H = A_.T(H, k)
            tot = tot + H
    return tot

def run(m, r, npts=4, seed=1):
    random.seed(seed + 100 * m + r)
    A_ = AHA(m)
    P = lhs_poly(A_, r)
    out = {}
    pts = []
    for _ in range(npts):
        X = [Fr(v, random.randint(1, 7)) for v in random.sample(range(2, 200), m)]
        X = list(dict.fromkeys(X)) if len(set(X)) == m else [Fr(v) for v in random.sample(range(2, 200), m)]
        s = Fr(random.randint(2, 9), random.randint(1, 5)); t = Fr(random.randint(2, 9), random.randint(1, 5))
        t = t if t != 1 else Fr(3, 2)
        pts.append((X, s, t))
    for conv in ('T', 'Tinv'):
        ratios = []
        for X, s, t in pts:
            L = evalp(A_, P, X, s, t); R = kernel(X, s, t, r, conv)
            ratios.append((L / R if R != 0 else None, t, s))
        out[conv] = ratios
    # operator identity check: sum_{W^J} T_w (pi^2 e_r) vs K_T
    F2 = A_.pi(A_.pi(A_.e(r)))
    S = parabolic_sym_op(A_, F2)
    opok = all(evalp(A_, S, X, s, t) == kernel(X, s, t, r, 'T') for X, s, t in pts)
    return out, opok, P, A_, pts

if __name__ == '__main__':
    for (m, r) in [(5, 2), (5, 3), (6, 2), (6, 3), (6, 4), (7, 5)]:
        out, opok, P, A_, pts = run(m, r)
        print(f'm={m} r={r}: operator sum_W^J T_w pi^2 e_r == K_T ? {opok}')
        for conv, rat in out.items():
            print(f'   conv={conv}: LHS/RHS at pts = ' + ', '.join(f'{x[0]} (t={x[1]},s={x[2]})' for x in rat))
        sys.stdout.flush()

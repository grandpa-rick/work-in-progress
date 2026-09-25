"""Day 206 step 4.
(a) Two-variable analogue of Lemma 2 (numerical, exact rationals):
    S_{n,p} := sum_{a != b} X_a^n X_b^p prod_{i in {a,b}, j notin {a,b}} (X_i - tX_j)/(X_i - X_j)
    Claim H:  (1-t)^2 S_{n,p} = (QJ(n,p) + QJ(p,n)) / (1+t),
    QJ(n,p) := [z^n w^p] Q(z)Q(w) F(w/z),  F(u) = (1-u)/(1-tu)   (Jing / Macdonald III (2.15))
            = sum_{k>=0} f_k q_{n+k} q_{p-k},  f_0 = 1, f_k = t^k - t^{k-1}.
    Equivalently (1-t)^2 R_(n,p) = QJ(n,p), where R_(n,p) = sum_{a!=b} X_a^n X_b^p a_ab prod_cross.
(b) Master formula: e_2 * e_r = t^{-1} e_2(Y).e_r = K(pi^2 e_r)
        = sum_{n,p} c_{np} QJ(n,p) / ((1+t)(1-t)^2),
    Phi(x,y) = xy [z^r] E(z)(1+sxz)(1+syz)/((1+xz)(1+yz)) = sum c_np x^n y^p, s = 1/q.
    Computed symbolically in free Lambda (e_1,e_2,... independent) and compared to the
    Day 191 conjecture W_r for r = 1..8.
"""
import itertools, random
from fractions import Fraction as Fr
import sympy as sp

# ---------- (a) numeric check of H ----------
def evals(vals, r):
    if r < 0 or r > len(vals): return Fr(0)
    tot = Fr(0)
    for c in itertools.combinations(vals, r):
        p = Fr(1)
        for v in c: p *= v
        tot += p
    return tot

def qn_vals(X, t, N):
    # Q(y) = E(-ty)/E(-y): q_n = sum_k (-t)^k e_k * h_{n-k}
    m = len(X); e = [evals(X, k) for k in range(N + 1)]
    h = [Fr(1)] + [Fr(0)] * N
    for n in range(1, N + 1):
        h[n] = sum((-1) ** (k + 1) * e[k] * h[n - k] for k in range(1, min(n, m) + 1))
    return [sum((-t) ** k * e[k] * h[n - k] for k in range(0, min(n, m) + 1)) for n in range(N + 1)]

def QJ_num(n, p, qn, t):
    tot = Fr(0)
    for k in range(0, p + 1):
        f = Fr(1) if k == 0 else t ** k - t ** (k - 1)
        tot += f * qn[n + k] * qn[p - k]
    return tot

def check_H(ms=(3, 4, 5, 6), maxdeg=5):
    random.seed(7); bad = 0; cnt = 0
    for m in ms:
        X = [Fr(v, random.randint(1, 5)) for v in random.sample(range(2, 99), m)]
        t = Fr(random.randint(2, 9), random.randint(1, 4)); t = t if t != 1 else Fr(5, 3)
        qn = qn_vals(X, t, 2 * maxdeg + 2)
        for n in range(1, maxdeg + 1):
            for p in range(1, maxdeg + 1):
                S = Fr(0); R = Fr(0)
                for a in range(m):
                    for b in range(m):
                        if a == b: continue
                        w = Fr(1)
                        for i in (a, b):
                            for j in range(m):
                                if j in (a, b): continue
                                w *= (X[i] - t * X[j]) / (X[i] - X[j])
                        S += X[a] ** n * X[b] ** p * w
                        R += X[a] ** n * X[b] ** p * w * (X[a] - t * X[b]) / (X[a] - X[b])
                lhsS = (1 - t) ** 2 * S
                rhsS = (QJ_num(n, p, qn, t) + QJ_num(p, n, qn, t)) / (1 + t)
                okR = (1 - t) ** 2 * R == QJ_num(n, p, qn, t)
                cnt += 1
                if lhsS != rhsS or not okR:
                    bad += 1; print('  H FAIL', m, n, p, lhsS == rhsS, okR)
    print(f'(a) Claim H [(1-t)^2 S_np = (QJ(n,p)+QJ(p,n))/(1+t)] and (1-t)^2 R_(n,p) = QJ(n,p):'
          f' {cnt - bad}/{cnt} OK  (m in {ms}, 1<=n,p<={maxdeg})')

# ---------- (b) symbolic master formula in free Lambda ----------
q, t, s = sp.symbols('q t s')
NMAX = 14
E = sp.symbols(f'e0:{NMAX + 1}')
def ee(k): return sp.Integer(0) if k < 0 else (sp.Integer(1) if k == 0 else E[k])
# h_n from e's, q_n = sum_k (-t)^k e_k h_{n-k}
H = [sp.Integer(1)]
for n in range(1, NMAX + 1):
    H.append(sp.expand(sum((-1) ** (k + 1) * ee(k) * H[n - k] for k in range(1, n + 1))))
QN = [sp.expand(sum((-t) ** k * ee(k) * H[n - k] for k in range(n + 1))) for n in range(NMAX + 1)]
def qq(n): return sp.Integer(0) if n < 0 else QN[n]
def QJ(n, p):
    return sp.expand(sum((sp.Integer(1) if k == 0 else t ** k - t ** (k - 1)) * qq(n + k) * qq(p - k)
                         for k in range(p + 1)))

def master(r):
    # Phi(x,y) = xy * sum_{i,j,u,v} e_{r-i-j-u-v} (-x)^i (-y)^j (s x)^u (s y)^v
    tot = sp.Integer(0)
    for i in range(r + 1):
        for j in range(r + 1 - i):
            for u in (0, 1):
                for v in (0, 1):
                    k = r - i - j - u - v
                    if k < 0: continue
                    c = ee(k) * (-1) ** (i + j) * s ** (u + v)
                    tot += c * QJ(1 + i + u, 1 + j + v)
    return sp.expand(tot / ((1 + t) * (1 - t) ** 2))

def qi(n): return sum(t ** i for i in range(n)) if n > 0 else 0
def W(r):
    return sp.expand(q ** -2 * ee(2) * ee(r) + (1 - 1 / q) / q * qi(r) * ee(1) * ee(r + 1)
                     + (1 - 1 / q) * sp.cancel(qi(r + 2) / qi(2)) * (qi(r + 1) - t * qi(r - 1) / q) * ee(r + 2))

if __name__ == '__main__':
    check_H()
    for r in range(1, 9):
        M = master(r).subs(s, 1 / q)
        d = sp.simplify(sp.together(M - W(r)))
        print(f'(b) r={r}: master formula - W_r (Day 191 conjecture) = {d}', flush=True)

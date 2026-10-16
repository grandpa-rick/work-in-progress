"""
Day 242 / Q473 referee check (Rick), independent of Clio's engine.

c -> infinity half-space functions of Garbali-de Gier-Mead-Wheeler, arXiv:2312.14348,
implemented from the PRINTED Theorem 5.1, eq. (5.1) (p.36-37), with
  Z_{L-n}(x_Kbar) -> prod_{j in Kbar} zeta(x_j),  zeta(x) = x(1-a x)/(x-a)   [Cor 4.7 / (5.13)]
  h(x) -> h_inf(x) = a(1-x^2)/(a-x)                                           [c -> oo of (1.2)/h]
  inner product prod_{j=1}^{nu_i - 1}  (as printed in (5.1); NOT the j-1 typo of (5.15)-(5.18))
nu = (nu_1 > nu_2 > ... > nu_n)  (nu_1 = right-most occupation, as in Lemma 5.2).

Arithmetic in F_p, p = 2^61-1, at a random (generic) point (q, a, y_1, ...).
"""
import random
from itertools import combinations, permutations

P = (1 << 61) - 1


def inv(x):
    x %= P
    if x == 0:
        raise ZeroDivisionError
    return pow(x, P - 2, P)


class Params:
    def __init__(self, seed, maxy=10, q=None, a=None, ys=None):
        r = random.Random(seed)
        self.q = r.randrange(2, P - 1) if q is None else q % P
        self.a = r.randrange(2, P - 1) if a is None else a % P
        self.y = {j: r.randrange(2, P - 1) for j in range(1, maxy + 1)}
        if ys:
            for j, v in ys.items():
                self.y[j] = v % P
        self.r = r

    def zeta(self, x):
        return x * (1 - self.a * x) % P * inv(x - self.a) % P

    def hinf(self, x):
        return self.a * (1 - x * x) % P * inv(self.a - x) % P

    def phi(self, m, x, shift=0):
        """phi_m(x); shift=1 plants an error: product runs to j <= m instead of m-1."""
        q, y = self.q, self.y
        v = (1 - q) * x % P * y[m] % P * inv(1 - q * x * y[m]) % P
        for j in range(1, m + shift):
            v = v * (1 - x * y[j]) % P * inv(1 - q * x * y[j]) % P
        return v


def G(par, kappa, xs, shift=0):
    """lim_{c->oo} G_kappa(x_1..x_L), kappa a decreasing tuple."""
    q = par.q
    kappa = tuple(sorted(kappa, reverse=True))
    d, L = len(kappa), len(xs)
    tot = 0
    for K in combinations(range(L), d):
        Kb = [j for j in range(L) if j not in K]
        pre = 1
        for j in Kb:
            pre = pre * par.zeta(xs[j]) % P
        for i in K:
            pre = pre * par.hinf(xs[i]) % P
            for j in Kb:
                xi, xj = xs[i], xs[j]
                pre = pre * (xj - q * xi) % P * inv(xj - xi) % P
                pre = pre * (1 - xi * xj) % P * inv(1 - q * xi * xj) % P
        for ii in range(d):
            for jj in range(ii + 1, d):
                xi, xj = xs[K[ii]], xs[K[jj]]
                pre = pre * (1 - xi * xj) % P * inv(1 - q * xi * xj) % P
        s = 0
        for sg in permutations(range(d)):
            t = 1
            for i in range(d):
                for j in range(i + 1, d):
                    xa, xb = xs[K[sg[j]]], xs[K[sg[i]]]
                    t = t * (xa - q * xb) % P * inv(xa - xb) % P
            for i in range(d):
                t = t * par.phi(kappa[i], xs[K[sg[i]]], shift) % P
            s = (s + t) % P
        tot = (tot + pre * s) % P
    return tot


def subsets(maxel, maxcard):
    out = []
    for d in range(maxcard + 1):
        for c in combinations(range(1, maxel + 1), d):
            out.append(tuple(sorted(c, reverse=True)))
    return out


def rank_modp(M):
    M = [row[:] for row in M]
    rk, cols = 0, len(M[0]) if M else 0
    for c in range(cols):
        piv = next((r for r in range(rk, len(M)) if M[r][c] % P), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        iv = inv(M[rk][c])
        M[rk] = [v * iv % P for v in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c] % P:
                f = M[r][c]
                M[r] = [(vr - f * vk) % P for vr, vk in zip(M[r], M[rk])]
        rk += 1
    return rk


def preceq(kp, k):
    """Clio Def 6: kp <= k iff |kp|<=|k| and kp_j <= k_j (j<=|kp|)."""
    return len(kp) <= len(k) and all(kp[j] <= k[j] for j in range(len(kp)))


def contained(kp, k):
    return set(kp) <= set(k)


def dominance_nocard(kp, k):
    """wrong order for control: entrywise on common length only, no cardinality clause."""
    m = min(len(kp), len(k))
    return all(kp[j] <= k[j] for j in range(m))

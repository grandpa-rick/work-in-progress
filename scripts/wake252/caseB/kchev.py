# Conjectural K-Chevalley for the (non-reduced) hyperplane section D = div(p_k) = sum m_nu X_nu on X_k in LG:
#   m_nu = 1 if k/nu is a diagonal box (a part 1 removed), else 2   (Chevalley multiplicity <omega_n, beta^vee>)
#   [O_D] = 1 - prod_nu (1 - m_nu x_nu + [m_nu=2] dX_nu),  products of Schubert classes by meet, dX_nu = reduced boundary of X_nu (Mobius)
# Check: H_k(j)-H_k(j-1) == sum_rho c_rho H_rho(j) at many rational j, all k with parts in [N].
import sys, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import gate as G
def m_of(k, nu): return 1 if len(nu) < len(k) else 2
def mul(A, B):
    C = defaultdict(int)
    for a, x in A.items():
        for b, y in B.items(): C[G.meet(a, b)] += x * y
    return {k: v for k, v in C.items() if v}
def kclass(k):
    P = {k: 1}   # k plays the role of 1 on X_k (meet with k is identity)
    for nu in G.covers_below(k):
        f = defaultdict(int); f[k] += 1; f[nu] -= m_of(k, nu)
        if m_of(k, nu) == 2:
            for rho, c in G.mobius_coeffs(nu).items(): f[rho] += c
        P = mul(P, {a: b for a, b in f.items() if b})
    out = defaultdict(int); out[k] += 1
    for rho, c in P.items(): out[rho] -= c
    return {r: c for r, c in out.items() if c}
def check(k, K):
    pts = [Fr(t, 2) for t in range(-4 * G.N - 4, 2 * sum(k) + 4)]
    return all(G.H(k, x) - G.H(k, x - 1) == sum(c * G.H(r, x) for r, c in K.items()) for x in pts)
if __name__ == '__main__':
    ok = tot = 0
    for L in range(1, G.N + 1):
        for c in itertools.combinations(range(1, G.N + 1), L):
            k = tuple(sorted(c, reverse=True)); K = kclass(k); t = check(k, K); ok += t; tot += 1
            print(k, t, {''.join(map(str, r)) or '0': v for r, v in K.items()}, flush=True)
    print(f'K-Chevalley guess: {ok}/{tot}')

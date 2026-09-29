"""Day 208: compute t^{-C(k,2)} e_k(Y).(e_a e_b) exactly, e-basis expansion.
Engine: day205 k3_fast_pipeline.AHA (conventions verbatim). e_k(Y) via the proved (A_k) per-tuple
identity Y_{b1}..Y_{bk}F = t^{C(k,2)} T_{w(b)} pi^k F, cross-checked against direct Y products (small m).
Output: dict partition -> coefficient as sympy expr in s,t (s = q^{-1}).
"""
import sys, time, itertools, pickle
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA, e_expansion
import sympy as sp
s_, t_ = sp.symbols('s t')

def sum_Tw(A, H, k):
    """sum over k-subsets d of [m] of T_{w(d)} H, w(d) = block_1 ... block_k, block_l = T_{d_l-1}..T_l
    (block_k acts first)."""
    m = A.m
    def rec(l, H, maxd):
        tot = 0 * A.s; V = H
        for d in range(l, maxd + 1):
            if d > l: V = A.T(V, d - 1)
            tot = tot + (V if l == 1 else rec(l - 1, V, d - 1))
        return tot
    return rec(k, H, m)

def ekY(A, F, k):
    """t^{-C(k,2)} e_k(Y) F"""
    H = F
    for _ in range(k): H = A.pi(H)
    return sum_Tw(A, H, k)

def ekY_direct(A, F, k):
    tot = 0 * A.s
    for b in itertools.combinations(range(1, A.m + 1), k):
        V = F
        for i in reversed(b): V = A.Y(V, i)
        tot = tot + V
    return tot

def to_st(A, c):
    m = A.m; ex = sp.Integer(0)
    for mon, co in c.to_dict().items():
        ex += sp.Integer(int(co)) * s_**mon[m] * t_**mon[m + 1]
    return ex

def compute(k, a, b, m=None):
    n = a + b + k
    if m is None: m = max(n, 1)
    A = AHA(m)
    F = A.e(a) * A.e(b)
    t0 = time.time()
    G = ekY(A, F, k)
    exp = e_expansion(A, G, n)
    return {lam: sp.expand(to_st(A, c)) for lam, c in exp.items()}, time.time() - t0

if __name__ == '__main__':
    # self-check: A_k-based ekY equals t^{-C(k,2)} * direct, small m
    for m in (3, 4, 5):
        A = AHA(m)
        for k in (1, 2, 3):
            if k > m: continue
            F = A.e(2) * A.e(1) + A.e(1) * A.e(1)
            d = ekY_direct(A, F, k) - A.t**(k*(k-1)//2) * ekY(A, F, k)
            assert d.is_zero(), (m, k)
    print('A_k engine == direct Y: OK (m<=5,k<=3)')

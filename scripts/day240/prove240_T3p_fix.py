"""T3' fixed: compose DFK Y_{i,1} symbolically at (q,theta), substitute (q,theta)->(s,theta^-1) ONCE at the end.
(prove240_check.py T3' substituted after each factor, which re-inverts earlier factors: script bug, not math.)"""
import sympy as sp
from prove240_check import symbasis, M_printed, ek_ops, z, mk, D_Yn, th, q, s, t
for N, deg in [(2, 3), (3, 3), (4, 2)]:
    x = mk(N); Fs = symbasis(x, deg)
    for k in range(1, N+1):
        b = sum(not z((q**(-k)*ek_ops(F, x, k, lambda g, i: D_Yn(g, x, i, 1))).subs({q: s, th: 1/th}, simultaneous=True)
                      - th**(k*(N-k))*M_printed(F, x, k, 1, s, 1/t)) for F in Fs)
        bc = sum(not z((q**(-k)*ek_ops(F, x, k, lambda g, i: D_Yn(g, x, i, 1))).subs({q: s}, simultaneous=True)
                      - th**(k*(N-k))*M_printed(F, x, k, 1, s, 1/t)) for F in Fs)
        print(f"N={N} k={k}: T3' DFK-ops = printed (1.5) at (s,theta^-1): bad={b}/{len(Fs)}   control (theta not inverted): bad={bc}", flush=True)

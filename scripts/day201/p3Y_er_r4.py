"""
Day 201: Compute p_3(Y) . e_4(X) at m = 7 to distinguish r-indep vs r-dep
coefficients in the DS-cone of (r, 1, 1, 1) = (4, 1, 1, 1).
"""
import sys, time
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day201')
from p3Y_er import compute_p3Y_er, partitions_of, dominance_ge
import sympy as sp

q, t = sp.symbols('q t')

r = 4
m = 7
print("=" * 72)
print(f"p_3(Y) . e_{r}(X) at m = {m}  (target degree n = {r + 3})")
print("=" * 72)
t0 = time.time()
exp = compute_p3Y_er(r, m)
dt = time.time() - t0
print(f"  compute time: {dt:.1f} s")

parts = partitions_of(r + 3)
pivot = tuple([r, 1, 1, 1])
print()
print("  e-basis expansion:")
for lam in parts:
    c = exp.get(lam, sp.Integer(0))
    cs = sp.factor(sp.simplify(c))
    print(f"    e_{lam}: {cs}")

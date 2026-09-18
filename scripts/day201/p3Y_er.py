"""
Day 201 Fallback (a): Compute p_3(Y) . e_r(X) for small r; extract e-basis expansion.

Meta-conjecture (Day 199 dream, level k = 3):
    p_3(Y) . e_r(X) has support size k+2 = 5 out of DS-cone of (r, 1, 1, 1),
    with k+1 = 4 r-independent coefficients plus one r-dependent
    tau^(3)_r . e_{r+3} where
        tau^(3)_r . q^5  =  A + B t^r + C t^{2r} + D t^{3r}   (Baxter-3).

Machinery copied from day198/p2Y_er.py (build_action). p_3(Y) = sum Y_i^3.
"""
import sympy as sp
from itertools import combinations
import time
import sys

# Add day198 to path so we can import
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day198')
from p2Y_er import (
    build_action, e_r_X, partitions_of, e_lambda_X,
    expand_symmetric_in_e_basis, dominance_ge,
)

q, t = sp.symbols('q t')


def compute_p3Y_er(r, m):
    """Return e-basis expansion dict of p_3(Y) . e_r(X)."""
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    erX = e_r_X(m, r)
    total = sp.Integer(0)
    for i in range(1, m+1):
        # Y_i^3 . e_r
        v = Y_apply(erX, i)
        v = Y_apply(v, i)
        v = Y_apply(v, i)
        total = sp.expand(total + v)
    n = r + 3
    return expand_symmetric_in_e_basis(total, m, n)


def report(r, m):
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
    print(f"  Pivot for DS-cone: {pivot}")
    print()
    print("  e-basis expansion:")
    for lam in parts:
        c = exp.get(lam, sp.Integer(0))
        cs = sp.factor(sp.simplify(c))
        is_zero = (sp.simplify(cs) == 0)
        ge_pivot = dominance_ge(lam, pivot)
        marker = ""
        if not is_zero and not ge_pivot:
            marker = "  <-- OUTSIDE DS-cone but NONZERO"
        if is_zero and ge_pivot:
            marker = "  (zero even though in DS-cone)"
        print(f"    e_{lam}  (>= pivot? {ge_pivot}, zero? {is_zero}):")
        print(f"      {cs}{marker}")
    return exp


if __name__ == "__main__":
    # r = 1 quick test: target degree 4, m = 4 minimum
    report(1, 4)
    # r = 2 at m = 5
    report(2, 5)
    # r = 3 at m = 6
    report(3, 6)

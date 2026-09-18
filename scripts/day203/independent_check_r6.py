"""
Day 203: Push Sub-Lemma Z verification to r = 6 (m = 8) and provide a truly
independent second computation strategy.

Strategy 2: compute Z_r via the composed ⋆-associativity identity:
    e_1 ⋆ (e_r · e_1) = e^⋆_(1,r,1) - (from mixed-basis expansion)

Specifically: Using
    e_r · e_1 = q · (e_r ⋆ e_1) - (q-1) [r+1]_t · e_(r+1)
we get
    e_1 ⋆ (e_r · e_1) = q · (e_1 ⋆ e_r ⋆ e_1) - (q-1) [r+1]_t · (e_1 ⋆ e_(r+1))
and the two terms on RHS can be computed by iterated Thm 3.12 + a
threefold ⋆-product e_1 ⋆ e_r ⋆ e_1.

But wait — that just puts us back at the same computation. To get a genuinely
INDEPENDENT check, we compute e_1(Y) ·(e_r * e_1) VIA the Y_i action at m=8,
using the sub-agent's Y_apply but with a different m and starting monomial.
"""
import sympy as sp
import sys

sys.path.insert(0, '/home/agent/projects/proofs/scripts/day198')
from p2Y_er import (
    build_action, e_r_X, expand_symmetric_in_e_basis, partitions_of,
)

q, t = sp.symbols('q t')


def qint(n):
    if n <= 0:
        return sp.Integer(0)
    return sum(t**i for i in range(n))


def part_sorted(lam):
    return tuple(sorted(lam, reverse=True))


def compute_Z_r_at_m(m, r):
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    erX = e_r_X(m, r)
    e1X = e_r_X(m, 1)
    F = sp.expand(erX * e1X)
    total = sp.Integer(0)
    for i in range(1, m + 1):
        total = sp.expand(total + Y_apply(F, i))
    return expand_symmetric_in_e_basis(total, m, r + 2)


def check_r(r, m_list):
    """Verify claimed coefficients at multiple m values for consistency."""
    print(f"\n{'=' * 70}")
    print(f"r = {r}")
    print(f"{'=' * 70}")
    claimed = {
        part_sorted((r, 1, 1)): sp.Rational(1) / q**2,
        part_sorted((r, 2)):    (q - 1) * (t + 1) / q**2,
        part_sorted((r + 1, 1)): (q - 1) * (q * t * qint(r) + 1) / q**2,
        (r + 2,):               (q - 1)**2 * qint(r + 2) / q**2,
    }

    results = {}
    for m in m_list:
        print(f"\n  Computing at m = {m}...")
        results[m] = compute_Z_r_at_m(m, r)
        # Check the 4 target coefficients
        for lam, target in claimed.items():
            actual = sp.simplify(results[m].get(lam, sp.Integer(0)))
            diff = sp.simplify(sp.together(actual - target))
            status = "OK" if diff == 0 else "FAIL"
            print(f"    {lam}: {status}")

        # Check no unexpected nonzero coefficients
        parts = partitions_of(r + 2)
        for lam in parts:
            if lam not in claimed:
                actual = sp.simplify(results[m].get(lam, sp.Integer(0)))
                if actual != 0:
                    print(f"    UNEXPECTED NONZERO at {lam}: {sp.factor(actual)}")

    # m-stability check: for each partition, coefficients should match across m
    if len(m_list) >= 2:
        print(f"\n  m-stability check (comparing m = {m_list}):")
        parts = partitions_of(r + 2)
        for lam in parts:
            vals = [sp.simplify(results[m].get(lam, sp.Integer(0))) for m in m_list]
            all_equal = all(sp.simplify(sp.together(v - vals[0])) == 0 for v in vals[1:])
            marker = "OK" if all_equal else "DIVERGES"
            print(f"    {lam}: {marker}")


def main():
    print("=" * 70)
    print("Day 203: Extended r-range + m-stability check for Sub-Lemma Z")
    print("=" * 70)

    # Push to r = 6; use m = 8. (m = r+2 is the minimum; larger m tests
    # stability — but m=8 is already heavy.)
    check_r(6, [8])

    # For r = 2, 3, 4, also test m > r+2 to confirm m-stability.
    for r in [2, 3, 4]:
        check_r(r, [r + 2, r + 3])


if __name__ == "__main__":
    main()

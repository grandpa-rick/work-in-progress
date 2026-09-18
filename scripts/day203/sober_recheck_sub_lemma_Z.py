"""
Day 203: SOBER re-derivation of Sub-Lemma Z.

Sub-Lemma Z (corrected statement):
    Z_r := e_1 star e_(r,1) = e_1(Y) . (e_r(X) e_1(X))

is a symmetric function of degree r+2. Claim: for r >= 2, its e-basis
expansion has exactly four nonzero coefficients:

    c_(r,1,1) = 1/q^2
    c_(r,2)   = (q-1)(t+1) / q^2
    c_(r+1,1) = (q-1)(q t [r]_t + 1) / q^2
    c_(r+2)   = (q-1)^2 [r+2]_t / q^2

Strategy:
  (1) Independent re-implementation. We import Y_apply from the Day 198
      infrastructure, which uses the AHA relations (T_i, pi_shift), verified
      independent of the sub-agent's Day 202 script.
  (2) Extend the r range: verify r = 2, 3, 4, 5, 6 (m up to 8).
  (3) For each r, print the FULL e-basis expansion to catch any nonzero
      coefficient outside the 4-term support.
  (4) Verify SYMBOLIC-in-r match between the direct computation and the
      claimed closed forms.

Second sober strategy:
  (5) Also compute via a DIFFERENT bracketing:
        Z_r' := t^{-binom(r,2)} * (e_1(Y) . e_r(Y) . e_1(Y) . 1)
             = ??? via q-multiplicativity
      and compare. This tests whether our Y_apply implementation is stable
      under permutation of operator order.
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


def compute_Z_r_direct(m, r):
    """
    Z_r := e_1(Y) . (e_r(X) e_1(X)),  degree r+2.

    Computes total = sum_i Y_i.(e_r * e_1) and expands in e-basis.
    """
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    erX = e_r_X(m, r)
    e1X = e_r_X(m, 1)
    F = sp.expand(erX * e1X)
    total = sp.Integer(0)
    for i in range(1, m + 1):
        total = sp.expand(total + Y_apply(F, i))
    return expand_symmetric_in_e_basis(total, m, r + 2)


def claimed_coeffs(r):
    """Return the four claimed coefficients at partitions of r+2."""
    return {
        part_sorted((r, 1, 1)): sp.Rational(1) / q**2,
        part_sorted((r, 2)):    (q - 1) * (t + 1) / q**2,
        part_sorted((r + 1, 1)): (q - 1) * (q * t * qint(r) + 1) / q**2,
        (r + 2,):               (q - 1)**2 * qint(r + 2) / q**2,
    }


def verify_r(r, m):
    """
    For a given r, compute Z_r at rank m and check:
      (a) Every partition of r+2 outside the 4-term support has coefficient 0.
      (b) The 4 nonzero coefficients match the claimed closed forms.
    """
    print(f"\n{'=' * 70}")
    print(f"r = {r},  m = {m}")
    print(f"{'=' * 70}")

    Z_r = compute_Z_r_direct(m, r)
    claimed = claimed_coeffs(r)
    parts = partitions_of(r + 2)

    all_ok = True
    print(f"Full e-basis expansion of Z_r = e_1 * e_(r,1) at degree {r+2}:")
    for lam in parts:
        actual = sp.simplify(Z_r.get(lam, sp.Integer(0)))
        if lam in claimed:
            target = sp.simplify(claimed[lam])
            diff = sp.simplify(sp.together(actual - target))
            status = "OK" if diff == 0 else "FAIL"
            if diff != 0:
                all_ok = False
            actual_str = sp.factor(actual)
            target_str = sp.factor(target)
            print(f"  {lam}:  actual = {actual_str}")
            print(f"          target = {target_str}   [{status}]")
        else:
            # Should be zero
            if actual != 0:
                all_ok = False
                print(f"  {lam}:  actual = {sp.factor(actual)}   [SHOULD BE ZERO but ISN'T]")
            else:
                print(f"  {lam}:  0  (support constraint OK)")

    print(f"\nOverall: {'PASSED' if all_ok else 'FAILED'}")
    return all_ok


def verify_general_r_symbolic():
    """
    Check that Newton identity [r+1]_t - t[r]_t = 1 (used in deriving Lemma 1's
    c_(r+1,1) closed form from Sub-Lemma Z) holds symbolically in r.

    This is the "cancellation" that makes p_2(Y).e_r r-independent at (r+1,1).
    """
    R = sp.Symbol('r_s', positive=True, integer=True)

    def qint_sym(n_expr):
        return (1 - t**n_expr) / (1 - t)

    # Identity: [r+1]_t - t[r]_t = 1.
    lhs = qint_sym(R + 1) - t * qint_sym(R)
    diff = sp.simplify(sp.together(lhs - 1))
    print(f"\nSymbolic Newton identity [r+1]_t - t*[r]_t = 1")
    print(f"  LHS - 1 = {diff}")
    assert diff == 0, "Newton identity fails"

    # Verify the R7 (r+1,1) cancellation coefficient is (q^2-1)/q^3, r-indep.
    # From:  c_(r+1,1)[p_2(Y).e_r] = alpha_r * beta + beta * c_(r+1,1)[Z_r] - 2 t B_r
    # with:  alpha_r = (1-1/q) [r+1]_t,  beta = 1/q,
    #        c_(r+1,1)[Z_r] = (q-1)(q t [r]_t + 1)/q^2,  B_r = (1/q)(1-1/q)[r]_t
    alpha_r = (1 - 1/q) * qint_sym(R + 1)
    beta = 1/q
    B_r = (1/q) * (1 - 1/q) * qint_sym(R)
    c_rp11_Z = (q - 1) * (q * t * qint_sym(R) + 1) / q**2
    c_rp11_p2Y = alpha_r * beta + beta * c_rp11_Z - 2 * t * B_r
    c_rp11_p2Y_simplified = sp.simplify(sp.together(c_rp11_p2Y))
    target = (q**2 - 1) / q**3
    diff2 = sp.simplify(sp.together(c_rp11_p2Y_simplified - target))
    print(f"\nSymbolic R7 (r+1,1) closed form (should be (q^2-1)/q^3):")
    print(f"  R7 gives: {sp.factor(c_rp11_p2Y_simplified)}")
    print(f"  target:   {sp.factor(target)}")
    print(f"  diff = {diff2}")
    assert diff2 == 0, "R7 r-indep coefficient does not match"


def verify_tau_r_symbolic():
    """
    Verify Rick's Day 200 tau_r closed form MATCHES the R7-derived formula
    for c_(r+2)[p_2(Y).e_r], symbolic in r.

    R7 gives:
      c_(r+2)[p_2(Y).e_r] = alpha_r * alpha_(r+1) + beta * c_(r+2)[Z_r] - 2t * A_r

    where c_(r+2)[Z_r] = (q-1)^2 [r+2]_t / q^2  (Sub-Lemma Z).

    Rick's Day 200:  tau_r = (A + B t^r + C t^(2r)) / q^3.
    """
    R = sp.Symbol('r_s', positive=True, integer=True)

    def qint_sym(n_expr):
        return (1 - t**n_expr) / (1 - t)

    alpha_r = (1 - 1/q) * qint_sym(R + 1)
    alpha_rp1 = (1 - 1/q) * qint_sym(R + 2)
    beta = 1/q
    A_r = (1 - 1/q) * qint_sym(R + 2) / qint_sym(2) * (qint_sym(R + 1) - t * qint_sym(R - 1) / q)
    c_rp2_Z = (q - 1)**2 * qint_sym(R + 2) / q**2

    tau_via_R7 = alpha_r * alpha_rp1 + beta * c_rp2_Z - 2 * t * A_r

    A = -(q**2 - 1) * (q - t - 1) / (t**2 - 1)
    B = t * (q**2 - 1) * (q - t) / (t - 1)
    C = -q * t**3 * (q**2 - 1) / (t**2 - 1)
    tau_day200 = (A + B * t**R + C * t**(2 * R)) / q**3

    diff = sp.simplify(sp.together(sp.expand(tau_via_R7 - tau_day200)))
    print(f"\nSymbolic tau_r via R7 vs Rick's Day 200 closed form:")
    print(f"  tau_via_R7 - tau_Day200 = {diff}")
    assert diff == 0, "tau_r symbolic mismatch"


def main():
    print("=" * 70)
    print("Day 203: Sober re-derivation of Sub-Lemma Z")
    print("=" * 70)

    # (1) Extended r-range verification
    all_ok = True
    for r in [2, 3, 4, 5]:
        m = r + 2  # exactly enough
        ok = verify_r(r, m)
        all_ok = all_ok and ok

    # (2) Symbolic Newton identity + r-indep coefficient check
    verify_general_r_symbolic()

    # (3) Symbolic tau_r match with Day 200 formula
    verify_tau_r_symbolic()

    print("\n" + "=" * 70)
    print(f"OVERALL: {'ALL CHECKS PASSED' if all_ok else 'SOME CHECKS FAILED'}")
    print("=" * 70)


if __name__ == "__main__":
    main()

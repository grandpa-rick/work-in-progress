"""
Day 202: Route R7 — Newton-cancellation attempt for Lemma 1 (p_2(Y)-Pieri).

Setup (Rick's Days 191, 198):
  - Hikita's Thm 3.12:  e_1 * e_r = (1 - q^-1) [r+1]_t e_{r+1} + q^-1 e_1 e_r.
  - Rick's e_2 * e_r (Day 191):  three-term closed form, r>=1.
  - Intertwiner property: e_a(Y) . F = t^binom(a,2) (e_a * F) for symmetric F.

Newton in Lambda(Y):
    p_2(Y) = e_1(Y)^2  -  2 e_2(Y).

Applied to e_r(X):
    p_2(Y) . e_r  =  e_1(Y).(e_1(Y).e_r)  -  2 e_2(Y).e_r
                  =  e_1(Y).(e_1 * e_r)   -  2 t (e_2 * e_r)          [since e_1*e_r symmetric]
                  =  e_1 * (e_1 * e_r)    -  2 t (e_2 * e_r).

So the "R7 Newton decomposition" is:
    p_2(Y) . e_r  =  e_1 * (e_1 * e_r)  -  2 t (e_2 * e_r).

The RHS uses only known closed forms (Thm 3.12 twice, Day 191 once).

Goal: verify this identity symbolically, and inspect where the cancellation happens.
"""
import sympy as sp

q, t = sp.symbols('q t')


def qint(n):
    """[n]_t = 1 + t + ... + t^{n-1}."""
    if n <= 0:
        return sp.Integer(0)
    return sum(t**i for i in range(n))


# =======================================================================
# Rick's / Hikita's closed forms in the ordinary e-basis.
# Represent each e-basis element by its partition tuple.
# A "Sym" is a dict {partition_tuple : coeff_in_Q(q,t)}.
# =======================================================================

def clean(D):
    """Drop zero entries, cancel each coefficient."""
    return {lam: sp.cancel(sp.together(c)) for lam, c in D.items()
            if sp.simplify(c) != 0}


def add(D1, D2):
    out = dict(D1)
    for k, v in D2.items():
        out[k] = out.get(k, 0) + v
    return clean(out)


def scale(D, s):
    return clean({k: v * s for k, v in D.items()})


def part_sorted(lam):
    """Return sorted partition tuple (weakly decreasing)."""
    return tuple(sorted(lam, reverse=True))


# ------------- Hikita Thm 3.12: e_1 * e_r ---------------

def e1_star_er(r):
    """
    e_1 * e_r = (1 - q^-1) [r+1]_t e_{r+1} + q^-1 e_{(r,1)}.
    Returned as {partition: coeff}.
    """
    out = {}
    coef_rp1 = (1 - 1/q) * qint(r + 1)
    coef_r1 = 1/q
    if r == 0:
        # e_1 * e_0 = e_1  (base)  -- but this shouldn't happen in our loop
        out[(1,)] = sp.Integer(1)
        return clean(out)
    out[(r+1,)] = coef_rp1
    out[part_sorted((r, 1))] = coef_r1
    return clean(out)


# ------------- Rick Day 191: e_2 * e_r ---------------

def e2_star_er(r):
    """
    e_2 * e_r = q^-2 e_{(r,2)}
              + q^-1 (1 - q^-1) [r]_t e_{(r+1,1)}
              + (1 - q^-1) [r+2]_t / [2]_t * ([r+1]_t - t [r-1]_t / q) e_{r+2}.
    Requires r >= 1.  For r=1 collapses via [0]_t = 0, [-1]_t = 0 conventions,
    but we may also cross-check against Thm 3.12 with commutativity.

    Note: for r=1 this must give e_2 * e_1 = e_1 * e_2 =
        (1-q^-1)[3]_t e_3 + q^-1 e_{(2,1)}.
    Let us verify.
    """
    out = {}
    # e_{(r,2)}
    p_r2 = part_sorted((r, 2))
    out[p_r2] = 1/q**2
    # e_{(r+1,1)}
    if r >= 1:
        c_rp11 = (1/q) * (1 - 1/q) * qint(r)
        out[part_sorted((r+1, 1))] = c_rp11
    # e_{r+2}
    if r >= 1:
        c_rp2 = (1 - 1/q) * qint(r + 2) / qint(2) \
                * (qint(r + 1) - t * qint(r - 1) / q)
        out[(r + 2,)] = c_rp2
    return clean(out)


# =======================================================================
# We need e_1 * (e_1 * e_r).  But e_1 * e_r is a linear combination of
# e_{r+1} and e_{(r,1)}.  We need e_1 * e_{r+1} AND e_1 * e_{(r,1)}.
#
# The second is NOT covered by Thm 3.12 directly!  It's e_1 * e_{lambda}
# for lambda = (r,1), a *composite* e-symbol.
#
# BUT: Rick's convention is that e_{lambda} = e_{lambda_1} * ... * e_{lambda_ell}
# is a Star-monomial (an element of Rick's Star-basis).  So e_{(r,1)}_STAR =
# e_r * e_1 (Star-product).  In the e-basis, however, e_{(r,1)} means the
# ordinary product e_r(X) e_1(X).  These are DIFFERENT.
#
# In Days 191/198, "e_{(r,1)}" in the RHS of Thm 3.12 refers to the ORDINARY
# product e_r(X) e_1(X) (symmetric function).  That's the e-basis element.
#
# So to compute e_1 * (e_1 * e_r) we need:
#     e_1 * e_{r+1}   (Thm 3.12 with r -> r+1)
#     e_1 * (e_r(X) e_1(X))   ← this is the depth-2 problem Rick's Day 191 flagged!
#
# By Day 198 Section 6, this is precisely the OBSTRUCTION.  Iterating Thm 3.12
# won't determine e_1 * (e_r e_1) without more input.  Rick's tautology is real.
#
# HOWEVER: we DO have an alternative.  e_1 * e_r * e_1 by ASSOCIATIVITY (Hikita
# 3.6) is the same three-way *-product regardless of parenthesization.  And
# Rick has C = e_1 * e_1 * e_2 computed at r=2 in the Day 198 verification.
#
# So we can COMPUTE e_1 * (e_1 * e_r) directly at r=2, 3, 4 via SymPy AHA,
# using the same p2Y_er.py framework.  Then check the identity
#
#     p_2(Y) . e_r  =  e_1 * (e_1 * e_r)  -  2t (e_2 * e_r)
#
# which is an *analytic* identity by Newton + intertwiner (both sides equal
# e_1(Y).(e_1(Y).e_r) - 2 e_2(Y).e_r; the RHS just re-expresses each atom via
# the intertwiner).
# =======================================================================

# Import Rick's AHA machinery
import sys
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day198')
from p2Y_er import build_action, e_r_X, expand_symmetric_in_e_basis, partitions_of


def compute_p2Y_er(m, r):
    """p_2(Y) . e_r(X) at rank m, returned as e-basis dict of degree r+2."""
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    erX = e_r_X(m, r)
    total = sp.Integer(0)
    for i in range(1, m+1):
        Yi_er = Y_apply(erX, i)
        Yi2_er = Y_apply(Yi_er, i)
        total = sp.expand(total + Yi2_er)
    return expand_symmetric_in_e_basis(total, m, r + 2)


def compute_e1_star_e1_star_er(m, r):
    """
    Compute C_r = e_1 * (e_1 * e_r) via direct AHA at rank m.

    Uses e_a(Y) . F = t^binom(a,2) (e_a * F) with a=1 (t^0=1).
    So e_1 * F = e_1(Y) . F for symmetric F.
    """
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    erX = e_r_X(m, r)
    # e_1(Y) . e_r
    e1Y_er = sp.Integer(0)
    for i in range(1, m+1):
        e1Y_er = sp.expand(e1Y_er + Y_apply(erX, i))
    # e_1(Y) . (e_1(Y) . e_r)
    e1Y_e1Y_er = sp.Integer(0)
    for i in range(1, m+1):
        e1Y_e1Y_er = sp.expand(e1Y_e1Y_er + Y_apply(e1Y_er, i))
    return expand_symmetric_in_e_basis(e1Y_e1Y_er, m, r + 2)


def compute_e2_star_er(m, r):
    """
    Compute e_2 * e_r via direct AHA: e_2 * F = t^-1 e_2(Y) . F.
    """
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    erX = e_r_X(m, r)
    # e_2(Y) . e_r  =  sum_{i<j} Y_i Y_j . e_r.  Iteratively Y_j then Y_i.
    total = sp.Integer(0)
    for i in range(1, m+1):
        for j in range(i+1, m+1):
            Yj_er = Y_apply(erX, j)
            YiYj_er = Y_apply(Yj_er, i)
            total = sp.expand(total + YiYj_er)
    # Now divide by t (the t^binom(2,2)=t^1 shift):  e_2 * e_r = t^-1 e_2(Y).e_r.
    total = sp.expand(total / t)
    return expand_symmetric_in_e_basis(total, m, r + 2)


def verify_R7(r):
    """Verify: p_2(Y).e_r == e_1*(e_1*e_r) - 2t (e_2*e_r) at rank m = r+2."""
    m = max(4, r + 2)
    print(f"\n=== Verify R7 at r={r} (rank m={m}) ===")
    p2Y_er = compute_p2Y_er(m, r)
    C_r = compute_e1_star_e1_star_er(m, r)
    e2se_r = compute_e2_star_er(m, r)

    # Compute RHS = C_r - 2t e2se_r  (as e-basis dict)
    rhs = {}
    for lam, c in C_r.items():
        rhs[lam] = rhs.get(lam, 0) + c
    for lam, c in e2se_r.items():
        rhs[lam] = rhs.get(lam, 0) - 2*t*c

    print(f"Partitions of {r+2}:")
    parts = partitions_of(r + 2)
    ok = True
    for lam in parts:
        lhs_c = sp.simplify(p2Y_er.get(lam, sp.Integer(0)))
        rhs_c = sp.simplify(rhs.get(lam, sp.Integer(0)))
        diff = sp.simplify(sp.expand(lhs_c - rhs_c))
        status = "OK" if diff == 0 else "FAIL"
        if diff != 0:
            ok = False
        print(f"  e_{lam}:  LHS={sp.factor(lhs_c)}   RHS={sp.factor(rhs_c)}   diff={diff}   [{status}]")
    return ok, p2Y_er, C_r, e2se_r


# =======================================================================
# The cancellation story: which coefficients are r-INDEP and how they arise.
# =======================================================================

def cancellation_story(rs):
    """
    Print, for each of the three r-INDEP positions (r,1,1), (r,2), (r+1,1),
    the contributions from C_r = e_1*(e_1*e_r) and from -2t (e_2*e_r).
    """
    print("\n\n" + "=" * 76)
    print("Cancellation story: what is r-dep / r-indep, term by term.")
    print("=" * 76)

    for r in rs:
        m = max(4, r + 2)
        p2Y_er = compute_p2Y_er(m, r)
        C_r = compute_e1_star_e1_star_er(m, r)
        e2se_r = compute_e2_star_er(m, r)

        # For r >= 2 the DS-relevant positions in degree r+2 are:
        #   (r+2), (r+1,1), (r,2), (r,1,1), and sub-dominant.
        positions = [(r+2,), part_sorted((r+1,1)), part_sorted((r,2)),
                     part_sorted((r,1,1))]
        # Only include ones with valid partitions of r+2
        positions = [p for p in positions if sum(p) == r + 2]

        print(f"\n-- r = {r}: contributions to each e-basis coefficient --")
        for lam in positions:
            c_lhs = sp.factor(sp.simplify(p2Y_er.get(lam, sp.Integer(0))))
            c_C = sp.factor(sp.simplify(C_r.get(lam, sp.Integer(0))))
            c_e2 = sp.factor(sp.simplify(e2se_r.get(lam, sp.Integer(0))))
            print(f"  e_{lam}:")
            print(f"     p_2(Y).e_r  = {c_lhs}")
            print(f"     e_1*e_1*e_r = {c_C}")
            print(f"     e_2*e_r     = {c_e2}")


def main():
    # Verify R7 for r = 2, 3, 4
    all_ok = True
    for r in [2, 3, 4]:
        ok, p2Y, C, e2 = verify_R7(r)
        if not ok:
            all_ok = False
    print()
    print("=" * 76)
    print(f"R7 Newton decomposition verified for r=2,3,4: {all_ok}")
    print("=" * 76)

    # Now the cancellation story
    cancellation_story([2, 3, 4])

    # Sub-leading r-INDEP check: show that
    #   c_(r,1,1) [p_2(Y).e_r]  is  r-independent (=1/q^3)
    #   c_(r,2)   [p_2(Y).e_r]  is  r-independent (= -(qt-q+t+1)/q^3)
    #   c_(r+1,1) [p_2(Y).e_r]  is  r-independent (= (q^2-1)/q^3)
    # BUT individually c_(r,1,1)[C_r] and c_(r,1,1)[2t e_2*e_r] are r-DEP.
    # Show explicitly.
    print("\n\n" + "=" * 76)
    print("Sub-leading r-INDEP check: C_r contribution at position (r,1,1)")
    print("=" * 76)
    for r in [2, 3, 4]:
        m = max(4, r + 2)
        C_r = compute_e1_star_e1_star_er(m, r)
        e2 = compute_e2_star_er(m, r)
        pos_r11 = part_sorted((r, 1, 1))
        c_C_r11 = sp.factor(sp.simplify(C_r.get(pos_r11, sp.Integer(0))))
        # c_(r,1,1) of 2t e_2*e_r should be 0 (support (r,1,1) not in SP interval)
        c_e2_r11 = sp.factor(sp.simplify(e2.get(pos_r11, sp.Integer(0))))
        print(f"  r={r}:  c_(r,1,1)[e_1*e_1*e_r] = {c_C_r11}   c_(r,1,1)[e_2*e_r] = {c_e2_r11}")


if __name__ == "__main__":
    main()

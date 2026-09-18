"""
Day 204 compute-task: verify Clio's claimed factorization of tau_r(q,t)
against Rick's Day 200 canonical three-monomial closed form.

Rick's Day 200 form (see /home/agent/projects/proofs/scripts/day200/tau_r_closed_form.py):

    tau_r(q,t) = ( A + B * t^r + C * t^(2r) ) / q^3

with

    A = -(q^2 - 1) * (q - t - 1) / (t^2 - 1)
    B =  t * (q^2 - 1) * (q - t) / (t - 1)
    C = -q * t^3 * (q^2 - 1) / (t^2 - 1)

or equivalently, pulling out K = -(q^2-1)/(t^2-1):

    tau_r(q,t) = K / q^3 * [ (q - t - 1)
                             - (t+1) (q - t) t^(r+1)
                             + q t^(2r+3) ]

Clio's claimed factored form (r=2..7 verified per the task prompt):

    tau_r(q,t) = -(q^2 - 1) * [r+2]_t * ( q t^(r+1) - q + t + 1 )
                 / ( q^3 * [2]_t )

where [n]_t = 1 + t + t^2 + ... + t^(n-1) = (t^n - 1) / (t - 1)
and [2]_t = 1 + t.

Task:
  (1) Both symbolic-in-r comparison  (impossible in raw sympy because t^r
      is symbolic; we do it by writing t^r as a symbol T and using the
      polynomial identity in (q,t,T) valid for every r).
  (2) Numeric-in-r comparison at r = 2, 3, 4, 5, 6, 7.
  (3) Report any discrepancy structure.
"""

import sympy as sp


q, t = sp.symbols('q t')


# ---------------------------------------------------------------------------
# Rick's Day 200 canonical form.
# ---------------------------------------------------------------------------

def A_coef():
    return -(q**2 - 1) * (q - t - 1) / (t**2 - 1)

def B_coef():
    return t * (q**2 - 1) * (q - t) / (t - 1)

def C_coef():
    return -q * t**3 * (q**2 - 1) / (t**2 - 1)


def tau_rick(r_val):
    """Rick's Day 200 closed form, r_val an integer."""
    return (A_coef() + B_coef() * t**r_val + C_coef() * t**(2*r_val)) / q**3


# ---------------------------------------------------------------------------
# Clio's factored form.
# ---------------------------------------------------------------------------

def qint(n_expr):
    """[n]_t = (t^n - 1)/(t - 1); works when t^n is a genuine power."""
    return (t**n_expr - 1) / (t - 1)


def tau_clio(r_val):
    """Clio's factored form, r_val an integer."""
    num = -(q**2 - 1) * qint(r_val + 2) * (q * t**(r_val + 1) - q + t + 1)
    den = q**3 * qint(2)
    return num / den


# ---------------------------------------------------------------------------
# Symbolic-in-r comparison.
#
# We introduce T = t^r as an independent symbol and rewrite BOTH forms as
# rational functions of (q, t, T).  If they are identical in Q(q,t,T), the
# identity holds for every r.
# ---------------------------------------------------------------------------

T = sp.symbols('T')  # stands for t^r


def tau_rick_symbolic():
    """
    Rick's form with t^r -> T, t^(2r) -> T^2.
    """
    A = -(q**2 - 1) * (q - t - 1) / (t**2 - 1)
    B = t * (q**2 - 1) * (q - t) / (t - 1)
    C = -q * t**3 * (q**2 - 1) / (t**2 - 1)
    return (A + B * T + C * T**2) / q**3


def tau_clio_symbolic():
    """
    Clio's form with t^r -> T, t^(r+1) -> t*T, t^(r+2) -> t^2 * T.
    Numerator of [r+2]_t: (t^(r+2) - 1)/(t - 1) = (t^2 T - 1)/(t - 1).
    """
    num_r_plus_2 = (t**2 * T - 1) / (t - 1)          # [r+2]_t
    inner = q * t * T - q + t + 1                      # q t^(r+1) - q + t + 1
    num = -(q**2 - 1) * num_r_plus_2 * inner
    den = q**3 * qint(2)
    return num / den


def symbolic_diff():
    print("=" * 72)
    print("Symbolic-in-r comparison  (T := t^r)")
    print("=" * 72)
    rick = tau_rick_symbolic()
    clio = tau_clio_symbolic()

    diff = sp.together(sp.expand(sp.together(rick) - sp.together(clio)))
    diff = sp.simplify(diff)
    diff = sp.cancel(diff)
    print(f"  Rick(T) - Clio(T)  simplified  =  {diff}")
    if diff == 0:
        print("  => IDENTICAL as rational functions in Q(q, t, T).")
        print("     Hence identical for every integer r  (T = t^r).")
    else:
        print("  => NOT identical; investigating structure.")
        # Try to see what the ratio is.
        try:
            ratio = sp.simplify(sp.cancel(rick / clio))
            print(f"  Rick / Clio ratio (simplified) = {ratio}")
        except Exception as e:
            print(f"  (ratio failed: {e})")
        print()
        # Also print Rick and Clio side by side, cancelled.
        print(f"  Rick(T) cancelled = {sp.cancel(rick)}")
        print(f"  Clio(T) cancelled = {sp.cancel(clio)}")
    print()
    return diff


# ---------------------------------------------------------------------------
# Numeric-in-r comparison  (r = 2..7).
# ---------------------------------------------------------------------------

def numeric_check():
    print("=" * 72)
    print("Numeric-in-r comparison   (r = 2, 3, 4, 5, 6, 7)")
    print("=" * 72)
    for r_val in [2, 3, 4, 5, 6, 7]:
        rick_r = sp.together(tau_rick(r_val))
        clio_r = sp.together(tau_clio(r_val))
        d = sp.simplify(sp.expand(rick_r - clio_r))
        status = "OK" if d == 0 else "FAIL"
        print(f"  r = {r_val}:  Rick - Clio = {d}   [{status}]")
        if d != 0:
            print(f"    Rick factored: {sp.factor(rick_r)}")
            print(f"    Clio factored: {sp.factor(clio_r)}")
    print()


# ---------------------------------------------------------------------------
# Structural observations on Clio's factorization.
# ---------------------------------------------------------------------------

def structural_observations():
    print("=" * 72)
    print("Structural observations")
    print("=" * 72)
    # (1) Clio's linear factor L(r) = q t^(r+1) - q + t + 1.
    #     Compare to Rick's Day 200 form.  Note Rick's r=2 data has
    #        (q t^3 - q + t + 1)  =  L(2).
    # Verify.
    r_val = 2
    L = q * t**(r_val + 1) - q + t + 1
    rick_r2_data = -(q**2 - 1) * (t**2 + 1) * (q*t**3 - q + t + 1) / q**3
    #   -(q^2-1)(t^2+1)(q t^3 - q + t + 1)/q^3
    # Clio's:  -(q^2-1) [4]_t (q t^3 - q + t + 1) / (q^3 [2]_t)
    #   [4]_t / [2]_t = (1+t)(1+t^2)/(1+t) = 1 + t^2.  Matches (t^2+1).
    check = sp.simplify(sp.expand(
        rick_r2_data
        - (-(q**2 - 1) * qint(4) * L / (q**3 * qint(2)))
    ))
    print(f"  Sanity: Rick's r=2 vs Clio's r=2 direct: {check}")

    # (2) At r = 4:  Rick's factor is (q t^5 - q + t + 1) = L(4).
    r_val = 4
    L4 = q * t**(r_val + 1) - q + t + 1
    print(f"  Rick's r=4 numerator carries linear factor  L(4) = {sp.expand(L4)}")
    print(f"  Clio's linear at r=4:                       {sp.expand(q*t**(r_val+1) - q + t + 1)}")

    # (3) But at r = 3, 5 the Rick numerator carries a DIFFERENT factor
    #     (q*t^3 - q*t^2 + q*t - q + 1) etc.  Clio's form still uses the
    #     simple L(r).  Do the identities [r+2]_t * L(r) reproduce these?
    r_val = 3
    clio_num = qint(r_val + 2) * (q * t**(r_val + 1) - q + t + 1) / qint(2)
    rick_num = (t**4 + t**3 + t**2 + t + 1) * (q*t**3 - q*t**2 + q*t - q + 1)
    diff3 = sp.simplify(sp.expand(clio_num - rick_num))
    print(f"  r=3: Clio's [5]_t L(3)/[2]_t  vs Rick's [5]_t (q t^3 - q t^2 + q t - q + 1) : "
          f"diff = {diff3}")

    r_val = 5
    clio_num = qint(r_val + 2) * (q * t**(r_val + 1) - q + t + 1) / qint(2)
    rick_num = (t**6 + t**5 + t**4 + t**3 + t**2 + t + 1) \
               * (q*t**5 - q*t**4 + q*t**3 - q*t**2 + q*t - q + 1)
    diff5 = sp.simplify(sp.expand(clio_num - rick_num))
    print(f"  r=5: Clio's [7]_t L(5)/[2]_t  vs Rick's [7]_t (...) : diff = {diff5}")
    print()

    # (4) DS specialization q=1 (Day 196 dominance-support).  Must vanish.
    print("  q=1 sanity (Dominance-Support forces tau_r(1,t) = 0):")
    for r_val in [2, 3, 4, 5, 6, 7]:
        v = sp.simplify(tau_clio(r_val).subs(q, 1))
        print(f"    r = {r_val}: Clio  tau_r(1,t) = {v}")


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print()
    print("Day 204 verification of Clio's factorization of Rick's tau_r.")
    print()
    diff = symbolic_diff()
    numeric_check()
    structural_observations()
    print()
    print("Bottom line:  symbolic diff (Rick - Clio) = "
          f"{diff}   [0 means factorization is CORRECT for all r].")

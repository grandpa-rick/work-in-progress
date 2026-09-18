"""
Day 204: Test Clio's Prediction 1 (Conjecture 10 at k=3).

Clio (2026-09-17 review of Rick's Day 200, §7.3):

    Prediction 1: q^5 * tau_r^(3) has (1 - t^{r+3}) = (1-t)[r+3]_t as a factor,
    and splits into three linear forms in u = t^r over Q(q,t).

Here tau_r^(3) := coefficient of e_{r+3} in p_3(Y).e_r(X)
(the "top" coefficient, analogous to k=1 with top at (r+1), k=2 with top at (r+2)).

Test methodology (over the rational function field Q(q)(t)):
  1. Take Rick's coefficient c_{(r+3)} of e_{(r+3)} from day201 logs.
  2. Multiply by q^5 to match Clio's normalization: q^5 * tau_r^(3).
  3. Compute polynomial division by [r+3]_t treating t as the variable,
     coefficients in Q(q). Report quot & rem.
  4. Also check: does q^5*tau_r^(3), regarded as a poly in t, vanish at
     each primitive d-th root of unity for every d | (r+3), d > 1?
     [r+3]_t = prod_{d | r+3, d > 1} Phi_d(t), so [r+3]_t divides iff
     the value vanishes at every such primitive root.
  5. If exact: try to write quotient as A0(t) + A1(t) u + A2(t) u^2 + A3(t) u^3
     with u = t^r and A_j(t) short, then factor as a cubic in u.
  6. Sanity: verify Clio's own k=2 formula divides by [r+2]_t.
"""
import sympy as sp
import cmath

q, t = sp.symbols('q t')
u = sp.Symbol('u')


def q_int(n):
    return sum(t**i for i in range(n))


# ------------------------------------------------------------------------
# TOP COEFFICIENTS tau_r^(3) at partition (r+3), verbatim from day201 logs
# /home/agent/projects/proofs/scripts/day201/p3Y_er{,_r4,_r5}.log
# ------------------------------------------------------------------------

tau3 = {}
tau3[1] = (q - 1)*(t + 1)*(t**2 + 1)*(q**2 + q + 1)*(q*t - q + 1)*(q**2*t**2 - q**2 + 1) / q**6
tau3[2] = ((q - 1)*(q**2 + q + 1)*(t**4 + t**3 + t**2 + t + 1) *
           (q**3*t**5 - q**3*t**4 - q**3*t + q**3 + q**2*t**3 - q**2 + q*t - q + 1)) / q**6
tau3[3] = ((q - 1)*(t + 1)*(q**2 + q + 1)*(t**2 - t + 1) *
           (q**3*t**9 - q**3*t**5 - q**3*t**4 + q**3 + q**2*t**6 + q**2*t**5
            + q**2*t**4 - q**2*t**2 - q**2*t - q**2 + q*t**3 - q + t**2 + t + 1)) / q**6
tau3[4] = ((q - 1)*(q**2 + q + 1)*(t**6 + t**5 + t**4 + t**3 + t**2 + t + 1) *
           (q**3*t**9 - q**3*t**8 + q**3*t**6 - q**3*t**5 - q**3*t**4
            + q**3*t**3 - q**3*t + q**3 + q**2*t**5 - q**2 + q*t - q + 1)) / q**6
tau3[5] = ((q - 1)*(t + 1)*(t**2 + 1)*(t**4 + 1)*(q**2 + q + 1) *
           (q**3*t**11 - q**3*t**10 + q**3*t**8 - q**3*t**7 - q**3*t**4
            + q**3*t**3 - q**3*t + q**3 + q**2*t**6 - q**2 + q*t - q + 1)) / q**6


# ------------------------------------------------------------------------
# Clio's k=2 top coefficient (§5.3):
#   q^3 tau_r = -(q^2-1)[r+2]_t(qt^{r+1}-q+t+1) / (q^3 [2]_t)  [box form]
# Equivalent:   = (q^2-1)(1-t^{r+2})(q-t-1-qt^{r+1})/(t^2-1)
# For sanity we use the second form since it makes the (1-t^{r+2}) factor
# manifest and lets us check rational-function divisibility by [r+2]_t.
# ------------------------------------------------------------------------

def clio_k2_top(r):
    return (q**2 - 1) * (1 - t**(r+2)) * (q - t - 1 - q * t**(r+1)) / (t**2 - 1)


# ------------------------------------------------------------------------
# Divisibility tests
# ------------------------------------------------------------------------

def poly_divide_by_qint(rational_expr, n):
    """
    Given a rational function `rational_expr` in t (with coeffs in Q(q)),
    check whether [n]_t = 1 + t + ... + t^{n-1} divides it.

    Strategy: write rational_expr as N(q,t)/D(q,t) with D not depending on t
    (or with D coprime to [n]_t). Then use sp.div on Poly(N, t).

    Returns (is_exact, quot_as_expr, rem_as_expr).
    """
    # Cancel to reduced form; separate polynomial in t.
    reduced = sp.cancel(rational_expr)
    N, D = sp.fraction(reduced)
    # If D still involves t, we may have a nontrivial denominator — check it
    # is coprime to [n]_t.
    qi = q_int(n)
    N_p = sp.Poly(sp.expand(N), t)
    qi_p = sp.Poly(qi, t)
    D_p = sp.Poly(sp.expand(D), t)
    # If D shares a factor with qi, the "divisibility" question is subtle.
    g = sp.gcd(D_p.as_expr(), qi_p.as_expr())
    if not (sp.Poly(g, t).degree() == 0):
        # denominator shares factor with [n]_t; adjust
        pass  # continue anyway — sp.div still meaningful in Q(q)[t]
    quot, rem = sp.div(N_p, qi_p, t)
    # Divide back by D
    return (rem.as_expr() == 0), sp.factor(quot.as_expr() / D), sp.factor(rem.as_expr() / D)


def vanishing_check(rational_expr, n):
    """
    Numerically check that rational_expr vanishes at all primitive d-th roots
    of unity for d | n, d > 1. Uses q = 1.7 as a generic value.
    Returns dict {d: |value|}.
    """
    results = {}
    qv = 1.7
    for d in sp.divisors(n):
        if d == 1:
            continue
        zeta = cmath.exp(2j * cmath.pi / d)
        # substitute t=zeta, q=qv numerically
        f = sp.lambdify((q, t), rational_expr, modules='numpy')
        try:
            val = complex(f(qv, zeta))
            results[d] = abs(val)
        except Exception as e:
            results[d] = f'error: {e}'
    return results


# ------------------------------------------------------------------------
# Structural extraction
# ------------------------------------------------------------------------

def extract_u_polynomial(rational_expr, r):
    """
    Given a rational function in t (over Q(q)), attempt to write it as
        A_0(t) + A_1(t) u + A_2(t) u^2 + A_3(t) u^3   with u = t^r
    where each A_j(t) is a polynomial of degree < r (a "residue").

    Returns dict {j: A_j(t)} and the u-polynomial as a sympy Expr in u,t,q.
    """
    reduced = sp.cancel(rational_expr)
    N, D = sp.fraction(reduced)
    poly = sp.Poly(sp.expand(N), t)
    coeffs = poly.all_coeffs()[::-1]  # index 0 = t^0
    A = {}
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        j, s = divmod(k, r)
        A.setdefault(j, sp.Integer(0))
        A[j] = sp.expand(A[j] + c * t**s)
    u_poly = sp.Add(*(sp.factor(A[j]) * u**j for j in sorted(A))) / D
    return A, u_poly


# ------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------

def main():
    print("=" * 74)
    print("Day 204: Test of Clio's Prediction 1 (Conjecture 10 at k=3)")
    print("=" * 74)

    # --- k=2 sanity check ---
    print("\n[SANITY k=2] q^3 * tau_r divisible by [r+2]_t (Clio §5.3):\n")
    for r in [2, 3, 4, 5]:
        expr = clio_k2_top(r)
        is_exact, quot, rem = poly_divide_by_qint(expr, r + 2)
        van = vanishing_check(expr, r + 2)
        status = "EXACT" if is_exact else "NON-EXACT"
        print(f"  r={r}: [{r+2}]_t divides q^3 tau_r? -> {status}")
        print(f"        vanishing at primitive d-th roots (d | {r+2}): {van}")
        if is_exact:
            print(f"        quot = {quot}")
    print()

    # --- k=3 prediction test ---
    print("=" * 74)
    print("[PREDICTION 1 k=3] q^5 * tau_r^(3) divisible by [r+3]_t ?")
    print("=" * 74)
    results_k3 = {}
    for r in [1, 2, 3, 4, 5]:
        expr = q**5 * tau3[r]
        n = r + 3
        print(f"\n--- r = {r},  target [n]_t with n = {n} ---")
        print(f"  q^5 * tau_r^(3) = {sp.factor(sp.together(expr))}")
        is_exact, quot, rem = poly_divide_by_qint(expr, n)
        van = vanishing_check(expr, n)
        print(f"  polynomial-division result: {'EXACT' if is_exact else 'NON-EXACT'}")
        print(f"  numerical vanishing at Phi_d for d | {n}: {van}")
        if is_exact:
            print(f"  quotient = {quot}")
        else:
            print(f"  remainder = {rem}")
        results_k3[r] = (is_exact, quot, rem)

    # --- Structural: u = t^r factorization ---
    print()
    print("=" * 74)
    print("[STRUCTURE] Extract A_j(t) so that quot = sum A_j(t) u^j (u = t^r)")
    print("=" * 74)
    for r in [2, 3, 4, 5]:
        is_exact, quot, rem = results_k3[r]
        if not is_exact:
            print(f"\n  r = {r}: division not exact; skipping u-decomposition")
            continue
        A, u_poly = extract_u_polynomial(quot, r)
        print(f"\n  r = {r}:")
        for j in sorted(A):
            print(f"    A_{j}(t) = {sp.factor(A[j])}")
        # Determine highest u-degree
        max_j = max(A) if A else 0
        print(f"    -> quot as poly in u (u=t^{r}) has degree {max_j}")
        print(f"       expression: {u_poly}")

    # --- Summary ---
    print()
    print("=" * 74)
    print("SUMMARY")
    print("=" * 74)
    print(f"{'r':>3} | {'n=r+3':>5} | {'[n]_t divides q^5 tau_r^(3)?':<32} | {'u-degree of quot':<20}")
    print("-" * 74)
    for r in [1, 2, 3, 4, 5]:
        is_exact, quot, rem = results_k3[r]
        if is_exact:
            A, _ = extract_u_polynomial(quot, r)
            max_j = max(A) if A else 0
            udeg = str(max_j)
        else:
            udeg = "N/A"
        print(f"{r:>3} | {r+3:>5} | {'YES (exact)' if is_exact else 'NO':<32} | {udeg:<20}")


if __name__ == "__main__":
    main()

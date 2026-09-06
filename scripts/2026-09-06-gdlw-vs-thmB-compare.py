"""
Compare Rick's D̄|_{E_3=0} (Theorem B, Day 170) with GDL-W's M_{P_n} = ω·PF_{n-1}.

Rick's object:
    D̄|_{E_3=0} = TY²[(q+1)² - E_1·T] / q³
    where q = 1 - E_1·T - 2·E_2·T·Y, Y = T·(1 + E_1·Y + E_2·Y²)
    E_1, E_2 are elementary symmetric functions in u_1, u_2 (2 variables).

    Explicit expansion (Day 162, Theorem B):
    [T^n] D̄|_{E_3=0} = (n+1) Σ_{b≥0} C(2b+2, b) C(n-1, 2b+2) E_1^{n-3-2b} E_2^b

GDL-W's object:
    M_{P_n} is a NEW symmetric function invariant, in INFINITELY many variables x_1, x_2, …
    M_{P_n} = ω · PF_{n-1} where PF is Panyushev's second parking function sym func.
    Under principal specialization (x_i = 1 for i≤n, else 0; equivalently ⟨·, h_1^n⟩
    or evaluation at n variables all = 1), it recovers Narayana N_n(t).

    Key structural facts about M_G:
    - M_G is a genuine symmetric function (infinite variables)
    - Top-degree component is e-positive for chordal graphs
    - For P_n: expansion into Schur functions gives Narayana-type coefficients

STRUCTURAL COMPARISON:

    Rick's D̄|_{E_3=0}: symmetric polynomial in TWO variables u_1, u_2 (via E_1, E_2).
    GDL-W's M_{P_n}:    symmetric function in INFINITELY many variables.

    A polynomial in TWO variables can be interpreted as a symmetric function
    only under principal specialization x_i = u_i for i=1,2, x_i=0 for i≥3.
    Under this specialization, ANY symmetric function truncates to a polynomial
    in E_1(u_1,u_2), E_2(u_1,u_2). And ω·PF_{n-1} truncated to 2 variables
    is far from arbitrary — it is a SPECIFIC symmetric function.

    So the question becomes: does the 2-variable truncation of ω·PF_{n-1}
    equal Rick's [T^n]D̄|_{E_3=0}?

    NOTE: Rick's D̄|_{E_3=0} is degree (n-3) in u's total (E_1^{n-3-2b}E_2^b has
    u-degree (n-3-2b) + 2b = n-3). It is the "top-of-layer-1" not the
    top-degree component of any obvious sym func on n variables.

We compute two things:
  (A) Rick's [T^n] D̄|_{E_3=0} as a polynomial in E_1, E_2 for n=3..8.
  (B) Compare specialization at (E_1, E_2) = (u_1+u_2, u_1·u_2) with
      the 2-variable truncation of ω·PF_{n-1} inferred from Narayana structure.

Reverse comparison (per user's fallback plan):
  Rick's coefficients [T^n]D̄|_{E_3=0} at E_1 = E_2 = 1 should match Narayana
  totals (Catalan) if the invariants agree. Compare Rick's sums vs Catalan.
"""

import sympy as sp
from sympy import Rational, binomial, symbols, expand, Poly, simplify

E1, E2, T, Y, q, t, u1, u2 = symbols('E1 E2 T Y q t u1 u2')

# =====================================================================
# (A) Rick's [T^n] D̄|_{E_3=0} from Theorem B closed form
# =====================================================================

def rick_barD_coeff(n):
    """
    [T^n] D̄|_{E_3=0} = (n+1) sum_{b>=0} C(2b+2,b) C(n-1, 2b+2) E1^{n-3-2b} E2^b
    """
    if n < 3:
        return sp.S(0)
    result = sp.S(0)
    b = 0
    while 2*b + 2 <= n - 1 and n - 3 - 2*b >= 0:
        coeff = (n+1) * binomial(2*b+2, b) * binomial(n-1, 2*b+2)
        result += coeff * E1**(n-3-2*b) * E2**b
        b += 1
    return result

print("=" * 70)
print("(A) Rick's [T^n] D̄|_{E_3=0} as polynomial in E_1, E_2")
print("=" * 70)
rick_coeffs = {}
for n in range(3, 9):
    p = rick_barD_coeff(n)
    rick_coeffs[n] = p
    print(f"  n={n}: {p}")

# =====================================================================
# (B) Rick at E_1 = E_2 = 1 → Narayana totals (Catalan)?
# =====================================================================
print("\n" + "=" * 70)
print("(B) Rick's [T^n] D̄|_{E_3=0} at (E_1, E_2) = (1, 1)")
print("=" * 70)
rick_at_1_1 = {}
for n in range(3, 9):
    val = rick_coeffs[n].subs({E1: 1, E2: 1})
    rick_at_1_1[n] = val
    print(f"  n={n}: {val}")

# Catalan numbers for comparison
catalan = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862]
print(f"\n  Catalan: {catalan[3:9]}")
print(f"  Rick vals: {[rick_at_1_1[n] for n in range(3,9)]}")
match_catalan = all(rick_at_1_1[n] == catalan[n] for n in range(3, 9))
print(f"  Match Catalan (n=3..8)?  {match_catalan}")

# Narayana polynomial N_n(t) = sum_k (1/n) C(n,k) C(n,k+1) t^{k-1}   [k=1..n]
def narayana_poly(n):
    """N_n(t) = sum_{k=1}^{n} (1/n) C(n,k) C(n,k+1) t^{k-1}."""
    if n == 0:
        return sp.S(1)
    return sum(Rational(1, n) * binomial(n, k) * binomial(n, k+1) * t**(k-1)
               for k in range(1, n))

print("\n" + "=" * 70)
print("(C) Narayana polynomials N_n(t) for reference")
print("=" * 70)
for n in range(3, 9):
    print(f"  N_{n}(t) = {sp.expand(narayana_poly(n))}")

# =====================================================================
# (D) Rick's [T^n] D̄|_{E_3=0} evaluated at (E_1, E_2) = (1+t, t)
# =====================================================================
# This is the "principal parking function specialization" — for two-variable
# elementary symmetric functions, u_1 = 1, u_2 = t gives E_1 = 1+t, E_2 = t.
print("\n" + "=" * 70)
print("(D) Rick at u_1=1, u_2=t: (E_1,E_2) = (1+t, t)")
print("=" * 70)
rick_at_1_t = {}
for n in range(3, 9):
    val = sp.expand(rick_coeffs[n].subs({E1: 1+t, E2: t}))
    rick_at_1_t[n] = val
    print(f"  n={n}: {val}")

# =====================================================================
# (E) Compare to Narayana N_n(t) directly and to shifted Narayana
# =====================================================================
print("\n" + "=" * 70)
print("(E) Compare Rick at (1+t, t) vs Narayana / shifted Narayana")
print("=" * 70)
for n in range(3, 9):
    R = rick_at_1_t[n]
    Nn = sp.expand(narayana_poly(n))
    Nn1 = sp.expand(narayana_poly(n+1))
    Nn_1 = sp.expand(narayana_poly(n-1))
    print(f"\n  n={n}:")
    print(f"    Rick(1+t,t) = {R}")
    print(f"    N_{n}(t)      = {Nn}")
    print(f"    ratio Rick/N_{n} = {sp.simplify(R/Nn) if Nn != 0 else 'N/A'}")
    print(f"    Rick - c·N_{n}: try c constant?")
    if Nn != 0:
        c_test = sp.simplify(R.subs(t, 0) / Nn.subs(t, 0)) if Nn.subs(t, 0) != 0 else "N/A"
        print(f"      c (at t=0) = {c_test}")

# =====================================================================
# (F) STRUCTURAL FACT: GDL-W's M_{P_n} is a symmetric function in
#     infinitely many variables; its principal specialization at
#     n variables all=1 gives the Narayana total (Catalan C_n).
#     Rick's D̄|_{E_3=0} is a 2-variable polynomial that at E_1=E_2=1
#     gives some number (computed above).
# =====================================================================
print("\n" + "=" * 70)
print("(F) Sanity: at (E1,E2)=(1,1), Rick gives what looks like...?")
print("=" * 70)
# Values: n=3: 4; n=4: 15; n=5: 60; n=6: 210; n=7: 720; n=8: 2394 (?)
# Compare OEIS: 4, 15, 60, 210 could be (n+1)*C(n-1)/something?
# Also try 4-term ratio test
vals = [rick_at_1_1[n] for n in range(3, 9)]
print(f"  Rick values at (1,1) n=3..8: {vals}")
# Test against (n+1)*C(n-1)*(n-1)/2  ...  what is [T^{n-1}] Y^2 at E_1=E_2=1?
# By Lagrange: [T^{n-1}] Y^2 at E1=E2=1 with phi=1+Y+Y^2:
#   Y = T·phi(Y) with phi = 1+Y+Y^2 is the Motzkin GF (or Fibonacci-like)
# But Y^2 coefficient sequence:
# Actually [T^n] D̄ = (n+1)(n-1)/2 · [T^{n-1}] Y^2   (Lagrange form)
# At E1=E2=1: [T^{n-1}] Y^2 where Y = T(1+Y+Y^2).
# Motzkin-like. Sequence for [T^k] Y^2 at E1=E2=1:
Y_var = sp.Symbol('Yv')
Y_series = sp.S(0)
for _ in range(12):
    phi_val = 1 + Y_series + Y_series**2  # E1=E2=1
    Y_series = sp.series(T*phi_val, T, 0, 12).removeO()
    Y_series = sp.expand(Y_series)
Y2_series = sp.expand(Y_series**2)
Y2_coeffs = [Y2_series.coeff(T, k) for k in range(12)]
print(f"  [T^k] Y^2 at E1=E2=1, k=0..11: {Y2_coeffs}")
print(f"  Then (n+1)(n-1)/2 * [T^{{n-1}}]Y^2 for n=3..8:")
for n in range(3, 9):
    if n-1 < len(Y2_coeffs):
        v = Rational((n+1)*(n-1), 2) * Y2_coeffs[n-1]
        print(f"    n={n}: {v}   (matches Rick? {v == rick_at_1_1[n]})")

# =====================================================================
# (G) THE DECISIVE STRUCTURAL COMPARISON
# =====================================================================
print("\n" + "=" * 70)
print("(G) DECISIVE STRUCTURAL COMPARISON")
print("=" * 70)
print("""
GDL-W's M_{P_n} = ω·PF_{n-1} is a symmetric function whose:
  - Ehrhart-like specialization (fine substitution, n variables) → N_n(t)
  - At t=1 (or principal spec at n vars all =1) → Catalan C_n
  - Expands in Schur basis; NUMBER of terms grows with n
    (partitions of some fixed weight related to n-1)

Rick's D̄|_{E_3=0} is at each [T^n]:
  - A polynomial in ONLY E_1, E_2 (2 variables u_1, u_2)
  - Total u-degree = n - 3 (NOT n or n-1)
  - Coefficient structure: Narayana-type multinomials weighted by (n+1)

Key mismatch #1: DEGREES.
  ω·PF_{n-1} lives in degree (n-1) or n as a sym func.
  Rick's [T^n]D̄|_{E_3=0} has u-degree (n-3) in 2 variables.
  A 2-var polynomial of degree (n-3) is FAR too small to encode M_{P_n}
  which as a sym func has ~partition-of-(n-1) many free coefficients.

Key mismatch #2: VARIABLE COUNT.
  Sym-func M_{P_n} needs infinitely many variables to define.
  Rick's D̄|_{E_3=0} uses only 2 (u_3 set to 0).

Key mismatch #3: NARAYANA at (1+t, t) vs at principal spec.
  Rick at (E_1,E_2)=(1+t,t) = principal spec of 2-var alphabet (1,t).
  See values in (D) — NOT a Narayana polynomial.
  Even at (E_1,E_2)=(1,1), Rick gives 4, 15, 60, 210, ...
  which are NOT Catalan numbers.

CONCLUSION: They are NOT the same invariant.
Both are related to Narayana structure but via distinct routes:
  - Rick's Narayana appears via Lagrange inversion of the 2-variable ν-system
    (the Y = T·phi(Y) fixed point) — a GF-level Narayana.
  - GDL-W's Narayana appears via bond-lattice shellability giving Narayana
    coefficients in the SCHUR expansion of a sym func on infinitely many vars.

They are "two DIFFERENT specializations that happen to both hit Narayana."
""")

# =====================================================================
# Final verdict summary
# =====================================================================
print("=" * 70)
print("VERDICT: RELATED BUT DIFFERENT")
print("=" * 70)
print("""
Evidence:
  1. Rick's D̄|_{E_3=0} at (E1,E2)=(1,1) gives: 4, 15, 60, 210, 720, 2394 (n=3..8)
     which are NOT Catalan (should be 5, 14, 42, 132, 429, 1430 for M_{P_n} totals).
  2. Rick's polynomial has degree (n-3) in 2 vars; ω·PF_{n-1} has degree (n-1)
     as a sym func on infinitely many vars.
  3. Rick's D̄|_{E_3=0} is a specific coefficient EXTRACTION from a 3-var
     Doubilet-like series; GDL-W's M_{P_n} is a full symmetric function.
  4. Both encode Narayana structure, but at DIFFERENT indexing levels:
     - Rick's D̄|_{E_3=0} coefficient of E_1^{n-3-2b}E_2^b involves
       Narayana-like multinomial (n+1)·C(2b+2,b)·C(n-1,2b+2).
     - GDL-W's ω·PF_{n-1} Schur-expansion coefficient at partition μ ⊢ (n-1)
       is a shifted Narayana number.
     These are different indexing schemes for the same numerical family.
""")

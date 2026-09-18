"""
Day 200: Closed form for tau_r(q,t), the r-dependent top coefficient of
p_2(Y) . e_r(X) expansion from Day 198.

Data (from Day 198 writeup / p2Y_er.py, p2Y_er_extended.py):

  tau_r(q,t) * q^3 =
    r=2: -(q^2-1)(t^2+1)(q t^3 - q + t + 1)
    r=3: -(q^2-1)(t^4+t^3+t^2+t+1)(q t^3 - q t^2 + q t - q + 1)
    r=4: -(q^2-1)(t^2-t+1)(t^2+t+1)(q t^5 - q + t + 1)
    r=5: -(q^2-1)(t^6+t^5+t^4+t^3+t^2+t+1)(q t^5 - q t^4 + q t^3 - q t^2 + q t - q + 1)

Clio's hint (email UID 276):
  (1-t) * (q [r]_t - t [r-2]_t) = (q-t) + t^(r-1)(1-qt)     ... "Baxterised"

Ansatz to test:
  tau_r(q,t) =? A(q,t) * [ (q-t) + t^(r-1)(1-qt) ] * B(q,t)
  tau_r(q,t) =? A [r]_t B + C t^(r-1) D   (linear-in-r-integer form)
  tau_r(q,t) =? A(q,t) + B(q,t) t^r + C(q,t) t^(2r)  (direct interpolation)
"""

import sympy as sp

q, t, r = sp.symbols('q t r')


def qint(n, x=t):
    """[n]_x = 1 + x + x^2 + ... + x^(n-1)."""
    if n == 0:
        return sp.Integer(0)
    return sum(x**i for i in range(n))


# ------------------------------------------------------------------
# tau_r data
# ------------------------------------------------------------------
def tau(r_val):
    """Return tau_r(q,t) for r in {2,3,4,5}."""
    if r_val == 2:
        return -(q**2 - 1) * (t**2 + 1) * (q*t**3 - q + t + 1) / q**3
    if r_val == 3:
        return -(q**2 - 1) * (t**4 + t**3 + t**2 + t + 1) \
               * (q*t**3 - q*t**2 + q*t - q + 1) / q**3
    if r_val == 4:
        return -(q**2 - 1) * (t**2 - t + 1) * (t**2 + t + 1) \
               * (q*t**5 - q + t + 1) / q**3
    if r_val == 5:
        return -(q**2 - 1) * (t**6 + t**5 + t**4 + t**3 + t**2 + t + 1) \
               * (q*t**5 - q*t**4 + q*t**3 - q*t**2 + q*t - q + 1) / q**3
    raise ValueError(r_val)


def report_tau_values():
    print("=" * 72)
    print("tau_r values (r = 2..5), factored:")
    print("=" * 72)
    for r_val in [2, 3, 4, 5]:
        tv = tau(r_val)
        print(f"  r = {r_val}: tau_r = {sp.factor(tv)}")


# ------------------------------------------------------------------
# Baxterised identity check
# ------------------------------------------------------------------
def check_baxterised_identity():
    print()
    print("=" * 72)
    print("Clio's Baxterised identity check:")
    print("  (1-t) (q [r]_t - t [r-2]_t) =?= (q-t) + t^(r-1)(1-qt)")
    print("=" * 72)
    for r_val in [2, 3, 4, 5, 6, 7]:
        lhs = (1 - t) * (q * qint(r_val) - t * qint(r_val - 2))
        rhs = (q - t) + t**(r_val - 1) * (1 - q*t)
        diff = sp.simplify(sp.expand(lhs - rhs))
        print(f"  r = {r_val}: LHS - RHS = {diff}")


# ------------------------------------------------------------------
# Divide by natural prefactor
# ------------------------------------------------------------------
def divide_by_prefactor():
    print()
    print("=" * 72)
    print("Divide tau_r by prefactor -(q^2-1)/q^3 = natural prefix.")
    print("Look at the residual U_r := tau_r * q^3 / [-(q^2-1)].")
    print("=" * 72)
    prefactor = -(q**2 - 1) / q**3
    for r_val in [2, 3, 4, 5]:
        U = sp.factor(sp.simplify(tau(r_val) / prefactor))
        print(f"  r = {r_val}: U_r = {U}")


# ------------------------------------------------------------------
# Structure of the "third factor" — the (q P(t) - Q(t) ...) piece
# ------------------------------------------------------------------
def third_factor():
    """
    From the data:
      r=2: third factor = q t^3 - q + t + 1
      r=3:                q t^3 - q t^2 + q t - q + 1
      r=4:                q t^5 - q + t + 1
      r=5:                q t^5 - q t^4 + q t^3 - q t^2 + q t - q + 1

    r=2 & r=4 have shape:  q t^{r+1} - q + t + 1  =  q(t^{r+1} - 1) + (t + 1)
    r=3 & r=5 have shape:  q * (alternating) + 1

    Let me try: F_r(q,t) := q t^{r+1} - q + t + 1 for r even
                F_r(q,t) := 1 + q * (t^{r+1} - ... )  for r odd (?)

    Alternative unified form: divide the second factor "S_r" from U_r/prefactor:
        r=2: S = (t^2+1), tail = (q t^3 - q + t + 1)
        r=3: S = [5]_t   , tail = (q t^3 - q t^2 + q t - q + 1)
        r=4: S = [3]_{t^2}, tail = (q t^5 - q + t + 1)
        r=5: S = [7]_t   , tail = (q t^5 - q t^4 + q t^3 - q t^2 + q t - q + 1)
    """
    print()
    print("=" * 72)
    print("Third-factor structure test")
    print("=" * 72)

    third_factors = {
        2: q*t**3 - q + t + 1,
        3: q*t**3 - q*t**2 + q*t - q + 1,
        4: q*t**5 - q + t + 1,
        5: q*t**5 - q*t**4 + q*t**3 - q*t**2 + q*t - q + 1,
    }

    for r_val, F in third_factors.items():
        # exponent seen: t^{r+1}
        Fs = sp.factor(F)
        print(f"  r = {r_val}: F_r = {sp.expand(F)}  factored = {Fs}")

    # Try: F_r = (t+1)(...) ?  or F_r as sum of q*P + Q form
    print()
    print("  Split F_r = q * P_r(t) + Q_r(t):")
    for r_val, F in third_factors.items():
        Fp = sp.Poly(F, q)
        P = Fp.coeff_monomial(q)      # coefficient of q
        Q = Fp.coeff_monomial(1)      # constant in q
        print(f"    r={r_val}: P_r(t) = {sp.factor(P)}   Q_r(t) = {sp.factor(Q)}")


# ------------------------------------------------------------------
# Direct Baxterised ansatz: does tau_r contain the factor
#   B_r := (q - t) + t^(r-1)(1 - q t)  ?
# ------------------------------------------------------------------
def baxterised_direct_test():
    print()
    print("=" * 72)
    print("Direct Baxterised test:")
    print("  Does B_r(q,t) := (q-t) + t^(r-1)(1-qt) divide tau_r(q,t) ?")
    print("=" * 72)

    for r_val in [2, 3, 4, 5]:
        B_r = (q - t) + t**(r_val - 1) * (1 - q*t)
        tau_val = tau(r_val)
        ratio = sp.together(tau_val / B_r)
        ratio_simplified = sp.cancel(ratio)
        # check if it's polynomial (numer/denom = ratio)
        num, den = sp.fraction(ratio_simplified)
        print(f"  r = {r_val}:")
        print(f"    B_r = {sp.expand(B_r)}")
        print(f"    tau_r / B_r = {sp.factor(ratio_simplified)}")
        # test whether B_r divides tau_r*q^3 as polynomials in q,t
        tau_num = sp.together(tau_val * q**3)
        tau_num = sp.cancel(tau_num).as_numer_denom()[0]
        try:
            quo, rem = sp.div(sp.Poly(tau_num, q, t), sp.Poly(B_r, q, t))
            print(f"    Polynomial division remainder = {rem.as_expr()}  (0 means B_r divides)")
        except Exception as e:
            print(f"    polynomial division error: {e}")


# ------------------------------------------------------------------
# Multiplicative ansatz using (1-t) (q [r]_t - t [r-2]_t)
# ------------------------------------------------------------------
def baxterised_qint_test():
    print()
    print("=" * 72)
    print("Baxterised q-integer form:")
    print("  B'_r := q [r]_t - t [r-2]_t  (dropping the (1-t))")
    print("  Test tau_r / B'_r for pattern.")
    print("=" * 72)

    for r_val in [2, 3, 4, 5]:
        B_prime = q * qint(r_val) - t * qint(r_val - 2)
        tau_val = tau(r_val)
        ratio = sp.cancel(sp.together(tau_val / B_prime))
        print(f"  r = {r_val}: B'_r = {sp.expand(B_prime)}")
        print(f"           tau_r / B'_r = {sp.factor(ratio)}")


# ------------------------------------------------------------------
# Ansatz:  tau_r(q,t) = A(q,t) + B(q,t) * t^r + C(q,t) * t^(2r)
# ------------------------------------------------------------------
def three_term_interp():
    print()
    print("=" * 72)
    print("Ansatz: tau_r * q^3 =?= A(q,t) + B(q,t) * t^r + C(q,t) * t^(2r).")
    print("Fit using r = 2,3,4 and predict r = 5. Then verify r = 5.")
    print("=" * 72)
    # We treat "unknowns" A,B,C in q, and known scalars t^r for r = 2,3,4.
    # Actually A,B,C are polynomials/rational functions in q,t independent of r.
    # 3 unknowns, 3 data points --> unique solve.
    A, B, C = sp.symbols('A B C')
    data = {r_val: sp.expand(tau(r_val) * q**3) for r_val in [2, 3, 4, 5]}

    eqns = []
    for r_val in [2, 3, 4]:
        # A + B t^r + C t^(2r) == data[r_val]
        eqns.append(A + B * t**r_val + C * t**(2*r_val) - data[r_val])

    sol = sp.solve(eqns, [A, B, C], dict=True)
    if not sol:
        print("  No solution.")
        return
    sol = sol[0]
    Aval = sp.factor(sp.simplify(sol[A]))
    Bval = sp.factor(sp.simplify(sol[B]))
    Cval = sp.factor(sp.simplify(sol[C]))
    print(f"  A(q,t) = {Aval}")
    print(f"  B(q,t) = {Bval}")
    print(f"  C(q,t) = {Cval}")

    # Check on r = 5
    pred5 = sp.simplify(Aval + Bval * t**5 + Cval * t**10)
    actual5 = sp.expand(tau(5) * q**3)
    diff5 = sp.simplify(sp.expand(pred5 - actual5))
    print()
    print(f"  Predicted (tau_5 * q^3) = {sp.factor(pred5)}")
    print(f"  Actual    (tau_5 * q^3) = {sp.factor(actual5)}")
    print(f"  Difference = {diff5}")
    print(f"  MATCH? {diff5 == 0}")


# ------------------------------------------------------------------
# Ansatz: tau_r*q^3 = A(q,t) + B(q,t) t^r  (2-term linear in t^r)
# ------------------------------------------------------------------
def two_term_interp():
    print()
    print("=" * 72)
    print("Ansatz: tau_r * q^3 =?= A(q,t) + B(q,t) t^r    (2-term).")
    print("=" * 72)
    A, B = sp.symbols('A B')
    data = {r_val: sp.expand(tau(r_val) * q**3) for r_val in [2, 3, 4, 5]}
    # Fit on r = 2, 3
    eqns = [A + B * t**2 - data[2],
            A + B * t**3 - data[3]]
    sol = sp.solve(eqns, [A, B], dict=True)
    if not sol:
        print("  No solution.")
        return
    sol = sol[0]
    Aval = sp.factor(sp.simplify(sol[A]))
    Bval = sp.factor(sp.simplify(sol[B]))
    print(f"  A(q,t) = {Aval}")
    print(f"  B(q,t) = {Bval}")

    for r_val in [4, 5]:
        pred = sp.simplify(Aval + Bval * t**r_val)
        actual = sp.expand(tau(r_val) * q**3)
        diff = sp.simplify(sp.expand(pred - actual))
        print(f"  r = {r_val}: diff = {sp.factor(diff)}  MATCH? {diff == 0}")


# ------------------------------------------------------------------
# Ansatz combining Baxterised and t^r pieces:
#   tau_r * q^3 = -(q^2-1) * [alpha(q,t) + beta(q,t)*t^(r-1) + gamma(q,t)*t^(2r)]
# ------------------------------------------------------------------
def refined_interp():
    print()
    print("=" * 72)
    print("Refined: tau_r*q^3 / (-(q^2-1)) = A(q,t) + B(q,t) t^r + C(q,t) t^(2r) ?")
    print("=" * 72)
    A, B, C = sp.symbols('A B C')
    pre = -(q**2 - 1)
    data = {r_val: sp.expand(sp.simplify(tau(r_val) * q**3 / pre)) for r_val in [2, 3, 4, 5]}
    eqns = [A + B * t**r_val + C * t**(2*r_val) - data[r_val] for r_val in [2, 3, 4]]
    sol = sp.solve(eqns, [A, B, C], dict=True)
    if not sol:
        print("  No solution.")
        return
    sol = sol[0]
    Aval = sp.factor(sp.simplify(sol[A]))
    Bval = sp.factor(sp.simplify(sol[B]))
    Cval = sp.factor(sp.simplify(sol[C]))
    print(f"  A(q,t) = {Aval}")
    print(f"  B(q,t) = {Bval}")
    print(f"  C(q,t) = {Cval}")

    pred5 = sp.simplify(Aval + Bval * t**5 + Cval * t**10)
    diff5 = sp.simplify(sp.expand(pred5 - data[5]))
    print(f"  r = 5 diff = {sp.factor(diff5)}  MATCH? {diff5 == 0}")


# ------------------------------------------------------------------
# Four-term ansatz: A + B t^r + C t^(r-1) + D t^(2r-1)
# (Motivated by Baxterised (q-t) + t^(r-1)(1-qt) structure.)
# ------------------------------------------------------------------
def four_term_interp():
    print()
    print("=" * 72)
    print("Ansatz: tau_r*q^3 = A + B t^r + C t^(2r) + D t^(r-1) (or similar 4-term).")
    print("Fit on r = 2,3,4,5.")
    print("=" * 72)
    A, B, C, D = sp.symbols('A B C D')
    data = {r_val: sp.expand(tau(r_val) * q**3) for r_val in [2, 3, 4, 5]}

    variants = [
        ("A + B t^r + C t^(2r) + D t^(2r-1)",
         lambda r_val: A + B*t**r_val + C*t**(2*r_val) + D*t**(2*r_val - 1)),
        ("A + B t^r + C t^(2r) + D t^(r+1)",
         lambda r_val: A + B*t**r_val + C*t**(2*r_val) + D*t**(r_val + 1)),
        ("A + B t^r + C t^(r-1) + D t^(2r+1)",
         lambda r_val: A + B*t**r_val + C*t**(r_val - 1) + D*t**(2*r_val + 1)),
    ]
    for name, ansatz in variants:
        print()
        print(f"  Ansatz: {name}")
        eqns = [ansatz(r_val) - data[r_val] for r_val in [2, 3, 4, 5]]
        sol = sp.solve(eqns, [A, B, C, D], dict=True)
        if not sol:
            print("    No solution (system inconsistent).")
            continue
        sol = sol[0]
        # need to check all variables solved
        Av = sp.factor(sp.simplify(sol.get(A, sp.nan)))
        Bv = sp.factor(sp.simplify(sol.get(B, sp.nan)))
        Cv = sp.factor(sp.simplify(sol.get(C, sp.nan)))
        Dv = sp.factor(sp.simplify(sol.get(D, sp.nan)))
        print(f"    A = {Av}")
        print(f"    B = {Bv}")
        print(f"    C = {Cv}")
        print(f"    D = {Dv}")


# ------------------------------------------------------------------
# The FULL Baxterised ansatz:
#   tau_r * q^3 = -(q^2-1) * G(q,t) * [(q-t) + t^(r-1)(1-qt)] * H_r(q,t)
# Try H_r = [r]_t, [r+2]_t, [r]_{t^2}, etc.
# ------------------------------------------------------------------
def full_baxterised_ansatz():
    print()
    print("=" * 72)
    print("Full Baxterised ansatz:")
    print("  tau_r * q^3 = -(q^2-1) * B_r * K_r(t)   where B_r=(q-t)+t^(r-1)(1-qt)")
    print("=" * 72)

    for r_val in [2, 3, 4, 5]:
        B_r = (q - t) + t**(r_val - 1) * (1 - q*t)
        pre = -(q**2 - 1)
        tau_val = tau(r_val)
        numer = sp.expand(tau_val * q**3)
        # tau_r * q^3 = pre * B_r * K   ==>   K = numer/(pre*B_r)
        K = sp.cancel(sp.together(numer / (pre * B_r)))
        num, den = sp.fraction(K)
        print(f"  r = {r_val}:")
        print(f"    B_r = {sp.expand(B_r)}")
        print(f"    K_r = tau_r * q^3 / (pre * B_r) = {sp.factor(K)}")


# ------------------------------------------------------------------
# Interpolation form using q-integers [k]_t
# tau_r * q^3 / [-(q^2-1)] = P(q,t) * [r+2]_t + Q(q,t) * t^r * (something)
# ------------------------------------------------------------------
def qint_split_test():
    """
    Observed: for odd r (3,5), the second factor is [r+2]_t.
              for even r (2,4), the second factor is [(r+2)/2]_{t^2} = (t^(r+2)-1)/(t^2-1).

    But  [r+2]_t = (t^(r+2)-1)/(t-1) while [(r+2)/2]_{t^2} = (t^(r+2)-1)/(t^2-1).

    Hmm.  Let's see:
       r=2 (even): [(r+2)/2]_{t^2} = [2]_{t^2} = 1 + t^2 = (t^4-1)/(t^2-1)
       r=4 (even): [3]_{t^2} = 1 + t^2 + t^4 = (t^6-1)/(t^2-1)
       r=3 (odd) : [5]_t = (t^5-1)/(t-1)
       r=5 (odd) : [7]_t = (t^7-1)/(t-1)

    The two shapes differ by a factor (t+1) between them:
       even: (t^(r+2)-1)/(t^2-1) = (t^(r+2)-1) / [(t-1)(t+1)]
       odd:  (t^(r+2)-1)/(t-1)

    Multiply even ones by (t+1):  even × (t+1) = [r+2]_t.
    So the "extra (t+1)" is soaked up somewhere else for odd r?

    Better idea: look at S_r * (third_factor) as one expression and see
    if the entire product (S_r * F_r) admits a cleaner unified form.
    """
    print()
    print("=" * 72)
    print("q-integer factorization audit: is the full product [S_r * F_r]")
    print("the natural atom?  Compute directly.")
    print("=" * 72)
    for r_val in [2, 3, 4, 5]:
        SF = sp.factor(sp.simplify(tau(r_val) * q**3 / (-(q**2 - 1))))
        SFe = sp.expand(SF)
        print(f"  r = {r_val}: S_r * F_r (expanded) = {SFe}")
        # Try to isolate (t^(r+2) - 1) as a factor:
        divisor = t**(r_val + 2) - 1
        try:
            quo, rem = sp.div(sp.Poly(SFe, q, t), sp.Poly(divisor, q, t))
            if rem.as_expr() == 0:
                print(f"    (t^(r+2)-1) divides!  Quotient = {sp.factor(quo.as_expr())}")
            else:
                print(f"    (t^(r+2)-1) does NOT divide as polynomial (rem = {sp.expand(rem.as_expr())})")
        except Exception as e:
            print(f"    poly div error: {e}")


# ------------------------------------------------------------------
# BEST FIT: after seeing S_r * F_r structure, try
#   tau_r * q^3 / (-(q^2-1)) = [ (t^(r+2)-1) * G(q,t) + t^(?) * H(q,t) ] / (t^2-1)?
# ------------------------------------------------------------------
def unified_ansatz():
    """
    Ansatz suggested by the data:
      tau_r * q^3 / (-(q^2-1)) = X(q,t) * (t^(r+2) - 1) / (t - 1)   ?
    Let's just check the numerator: does (t^(r+2)-1) divide? (see qint_split_test above).
    If not, try:
      tau_r * q^3 = -(q^2-1)/(t-1) * [ alpha(q,t) * t^(2r+something) + ... ]
    """
    print()
    print("=" * 72)
    print("Unified interpolation: solve for A(q,t), B(q,t), C(q,t), D(q,t) satisfying")
    print("  tau_r * q^3 * (t-1) / (-(q^2-1)) =?= A + B t^r + C t^(2r) + D t^(2r+1)")
    print("=" * 72)

    A, B, C, D = sp.symbols('A B C D')
    pre = -(q**2 - 1) / (t - 1)
    data = {r_val: sp.expand(sp.simplify(tau(r_val) * q**3 / pre)) for r_val in [2, 3, 4, 5]}
    eqns = [A + B*t**r_val + C*t**(2*r_val) + D*t**(2*r_val+1) - data[r_val]
            for r_val in [2, 3, 4, 5]]
    sol = sp.solve(eqns, [A, B, C, D], dict=True)
    if not sol:
        print("  No solution.")
        return
    sol = sol[0]
    for name, sym in [('A', A), ('B', B), ('C', C), ('D', D)]:
        val = sp.factor(sp.simplify(sol.get(sym, sp.nan)))
        print(f"  {name} = {val}")


# ------------------------------------------------------------------
# Even simpler: is tau_r/(-(q^2-1)/q^3) = polynomial in q, t^r, and t^(2r)?
# ------------------------------------------------------------------
def look_at_expanded_form():
    print()
    print("=" * 72)
    print("Expanded form of U_r := tau_r * q^3 / (-(q^2-1)):")
    print("=" * 72)
    pre = -(q**2 - 1)
    for r_val in [2, 3, 4, 5]:
        U = sp.expand(sp.simplify(tau(r_val) * q**3 / pre))
        # Extract q-degrees
        Up = sp.Poly(U, q)
        deg_q = Up.degree()
        print(f"  r = {r_val}: q-degree = {deg_q}")
        for d in range(deg_q + 1):
            coef = sp.factor(Up.coeff_monomial(q**d))
            print(f"    q^{d} coeff = {coef}")


# ------------------------------------------------------------------
# Direct linear ansatz over t-monomials:
# split U_r into t-degrees, see which exponents actually appear as r varies.
# ------------------------------------------------------------------
def t_degree_scan():
    print()
    print("=" * 72)
    print("t-degree pattern of tau_r * q^3 / (-(q^2-1)):")
    print("=" * 72)
    pre = -(q**2 - 1)
    for r_val in [2, 3, 4, 5]:
        U = sp.expand(sp.simplify(tau(r_val) * q**3 / pre))
        # Get all (q_deg, t_deg) monomials with coefficients
        # coefficients are integers/rationals.
        # We treat U as polynomial in q, t.
        Uc = sp.Poly(U, q, t)
        monoms = Uc.monoms()
        coeffs = Uc.coeffs()
        # Group by t-degree
        by_t = {}
        for m, c in zip(monoms, coeffs):
            qd, td = m
            by_t.setdefault(td, []).append((qd, c))
        print(f"  r = {r_val}: t-degrees present = {sorted(by_t.keys())}")
        for td in sorted(by_t.keys()):
            desc = ", ".join(f"{c}*q^{qd}" for qd, c in by_t[td])
            print(f"    t^{td}: {desc}")


# ------------------------------------------------------------------
# Try the "even/odd unified" ansatz:
#   U_r := tau_r * q^3 / (-(q^2-1)) = c(q,t) * [r+2]_t + d(q,t) * t^r * [something]
# ------------------------------------------------------------------
def qint_r2_ansatz():
    print()
    print("=" * 72)
    print("Ansatz: U_r =?= alpha(q,t) [r+2]_t + beta(q,t) t^r [r]_t + ...")
    print("Try:  U_r = [ f1(q,t) + f2(q,t) t^r + f3(q,t) t^(2r+2) ] / (t-1)")
    print("=" * 72)
    A, B, C = sp.symbols('A B C')
    pre = -(q**2 - 1)
    data = {r_val: sp.expand(sp.simplify(tau(r_val) * q**3 / pre * (t - 1)))
            for r_val in [2, 3, 4, 5]}
    for r_val in [2, 3, 4, 5]:
        print(f"    (t-1)*U_{r_val} expanded = {data[r_val]}")

    # Try ansatz: (t-1) U_r = A + B t^r + C t^(2r+2)
    eqns = [A + B*t**r_val + C*t**(2*r_val + 2) - data[r_val] for r_val in [2, 3, 4]]
    sol = sp.solve(eqns, [A, B, C], dict=True)
    if sol:
        sol = sol[0]
        Aval = sp.factor(sp.simplify(sol[A]))
        Bval = sp.factor(sp.simplify(sol[B]))
        Cval = sp.factor(sp.simplify(sol[C]))
        print()
        print(f"  Fit on r = 2,3,4:")
        print(f"    A = {Aval}")
        print(f"    B = {Bval}")
        print(f"    C = {Cval}")
        # check r=5
        pred5 = sp.simplify(Aval + Bval*t**5 + Cval*t**12)
        diff5 = sp.simplify(sp.expand(pred5 - data[5]))
        print(f"    r = 5 diff = {sp.factor(diff5)}  MATCH? {diff5 == 0}")


# ------------------------------------------------------------------
# Cleanest form to try:  (t-1) * U_r = P + Q * t^(2r+2)  (2-term!)
# Motivated by looking at t-degree scan structure.
# ------------------------------------------------------------------
def cleanest_two_term():
    print()
    print("=" * 72)
    print("Try: (t-1) * U_r = P(q,t) + Q(q,t) t^(2r+2)")
    print("=" * 72)
    A, B = sp.symbols('A B')
    pre = -(q**2 - 1)
    data = {r_val: sp.expand(sp.simplify(tau(r_val) * q**3 / pre * (t - 1)))
            for r_val in [2, 3, 4, 5]}
    eqns = [A + B*t**(2*r_val + 2) - data[r_val] for r_val in [2, 3]]
    sol = sp.solve(eqns, [A, B], dict=True)
    if not sol:
        print("  no solution")
        return
    sol = sol[0]
    Aval = sp.factor(sp.simplify(sol[A]))
    Bval = sp.factor(sp.simplify(sol[B]))
    print(f"  A = {Aval}")
    print(f"  B = {Bval}")
    for r_val in [4, 5]:
        pred = sp.simplify(Aval + Bval*t**(2*r_val + 2))
        diff = sp.simplify(sp.expand(pred - data[r_val]))
        print(f"  r = {r_val}: MATCH? {diff == 0}")


# ------------------------------------------------------------------
# Final serious attempt: A + B t^r + C t^(r+1) + D t^(2r+2)
# ------------------------------------------------------------------
def four_term_ansatz_v2():
    print()
    print("=" * 72)
    print("Ansatz: (t-1)*U_r = A + B t^r + C t^(r+1) + D t^(2r+2)")
    print("=" * 72)
    A, B, C, D = sp.symbols('A B C D')
    pre = -(q**2 - 1)
    data = {r_val: sp.expand(sp.simplify(tau(r_val) * q**3 / pre * (t - 1)))
            for r_val in [2, 3, 4, 5]}
    eqns = [A + B*t**r_val + C*t**(r_val+1) + D*t**(2*r_val + 2) - data[r_val]
            for r_val in [2, 3, 4, 5]]
    sol = sp.solve(eqns, [A, B, C, D], dict=True)
    if not sol:
        print("  system inconsistent, no fit")
        return
    sol = sol[0]
    for name, sym in [('A', A), ('B', B), ('C', C), ('D', D)]:
        val = sp.factor(sp.simplify(sol.get(sym, sp.nan)))
        print(f"  {name} = {val}")


# ------------------------------------------------------------------
# main
# ------------------------------------------------------------------
if __name__ == "__main__":
    report_tau_values()
    check_baxterised_identity()
    divide_by_prefactor()
    third_factor()
    baxterised_direct_test()
    baxterised_qint_test()
    two_term_interp()
    three_term_interp()
    refined_interp()
    four_term_interp()
    full_baxterised_ansatz()
    qint_split_test()
    unified_ansatz()
    look_at_expanded_form()
    t_degree_scan()
    qint_r2_ansatz()
    cleanest_two_term()
    four_term_ansatz_v2()

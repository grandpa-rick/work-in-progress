"""
Day 202: Independent derivation of Rick's Day 200 tau_r closed form
via R7 + Sub-Lemma Z at position (r+2).

We use:
    c_(r+2)[p_2(Y).e_r]  =  alpha_r alpha_{r+1}  +  beta * c_(r+2)[Z_r]  -  2 t A_r

with
    alpha_r = (1 - 1/q) [r+1]_t
    alpha_{r+1} = (1 - 1/q) [r+2]_t
    beta = 1/q
    A_r = (1 - 1/q) [r+2]_t / [2]_t * ([r+1]_t - t [r-1]_t / q)
    c_(r+2)[Z_r] = (q-1)^2 [r+2]_t / q^2      (Sub-Lemma Z at position 4)

Compare against Rick's Day 200 closed form:
    tau_r = (A + B t^r + C t^{2r}) / q^3
    A = -(q^2-1)(q - t - 1) / (t^2 - 1)
    B =  t (q^2-1)(q - t) / (t - 1)
    C = -q t^3 (q^2-1) / (t^2 - 1)
"""
import sympy as sp

q, t = sp.symbols('q t')


def qint(n):
    if n <= 0:
        return sp.Integer(0)
    return sum(t**i for i in range(n))


def tau_via_R7(r):
    """Compute c_(r+2)[p_2(Y).e_r] via R7 + Sub-Lemma Z position 4."""
    alpha_r = (1 - 1/q) * qint(r + 1)
    alpha_rp1 = (1 - 1/q) * qint(r + 2)
    beta = 1/q
    A_r = (1 - 1/q) * qint(r + 2) / qint(2) * (qint(r + 1) - t * qint(r - 1) / q)
    c_rp2_Zr = (q - 1)**2 * qint(r + 2) / q**2

    val = alpha_r * alpha_rp1 + beta * c_rp2_Zr - 2*t*A_r
    return sp.simplify(sp.expand(val))


def tau_day200(r):
    """Rick's Day 200 closed form."""
    A = -(q**2 - 1) * (q - t - 1) / (t**2 - 1)
    B = t * (q**2 - 1) * (q - t) / (t - 1)
    C = -q * t**3 * (q**2 - 1) / (t**2 - 1)
    return sp.simplify((A + B * t**r + C * t**(2*r)) / q**3)


def main():
    print("=" * 72)
    print("Independent derivation of Rick's tau_r via R7 + Sub-Lemma Z(r+2)")
    print("=" * 72)
    for r in [2, 3, 4, 5, 6]:
        via_R7 = tau_via_R7(r)
        day200 = tau_day200(r)
        diff = sp.simplify(sp.expand(via_R7 - day200))
        status = "OK" if diff == 0 else "FAIL"
        print(f"\n  r = {r}:")
        print(f"    R7:      {sp.factor(via_R7)}")
        print(f"    Day 200: {sp.factor(day200)}")
        print(f"    diff:    {diff}   [{status}]")

    # General r check: verify symbolically that
    #   [(1-1/q)^2 [r+1]_t [r+2]_t] + [(q-1)^2 [r+2]_t / q^3]
    #     - [2t (1-1/q) [r+2]_t / [2]_t * ([r+1]_t - t [r-1]_t / q)]
    #   equals Rick's (A + B t^r + C t^{2r}) / q^3.
    #
    # We do this by casting r as a symbol and using [n]_t = (1 - t^n)/(1 - t).
    print("\n\n" + "=" * 72)
    print("General r: analytic equality via [n]_t = (1 - t^n)/(1 - t).")
    print("=" * 72)
    R = sp.Symbol('r_sym', integer=True, positive=True)

    def qint_sym(n_expr):
        """(1 - t^n)/(1 - t)."""
        return (1 - t**n_expr) / (1 - t)

    alpha_r_sym = (1 - 1/q) * qint_sym(R + 1)
    alpha_rp1_sym = (1 - 1/q) * qint_sym(R + 2)
    beta = 1/q
    A_r_sym = (1 - 1/q) * qint_sym(R + 2) / qint_sym(2) * (qint_sym(R + 1) - t * qint_sym(R - 1) / q)
    c_rp2_Zr_sym = (q - 1)**2 * qint_sym(R + 2) / q**2

    tau_sym = alpha_r_sym * alpha_rp1_sym + beta * c_rp2_Zr_sym - 2*t*A_r_sym

    A_coef = -(q**2 - 1) * (q - t - 1) / (t**2 - 1)
    B_coef = t * (q**2 - 1) * (q - t) / (t - 1)
    C_coef = -q * t**3 * (q**2 - 1) / (t**2 - 1)
    tau_day200_sym = (A_coef + B_coef * t**R + C_coef * t**(2*R)) / q**3

    diff = sp.simplify(sp.expand(sp.together(tau_sym - tau_day200_sym)))
    print(f"  tau_via_R7 - tau_Day200 (symbolic in r) = {diff}")
    if diff == 0:
        print("  SYMBOLIC IDENTITY VERIFIED for all r.")
    else:
        print(f"  simplify further ...")
        diff2 = sp.simplify(diff.rewrite(sp.Piecewise))
        print(f"    -> {diff2}")


if __name__ == "__main__":
    main()

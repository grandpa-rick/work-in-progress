"""
Day 200: CLEAN closed form for tau_r(q,t).

Discovery (this session):

    tau_r(q,t) * q^3  =  A(q,t) + B(q,t) t^r + C(q,t) t^(2r)

with r-independent coefficients

    A(q,t) = -(q^2-1)(q - t - 1) / (t^2 - 1)      = -(q-1)(q+1)(q-t-1)/((t-1)(t+1))
    B(q,t) =  t (q^2-1)(q - t) / (t - 1)          =  t(q-1)(q+1)(q-t)/(t-1)
    C(q,t) = -q t^3 (q^2-1) / (t^2 - 1)           = -q t^3 (q-1)(q+1)/((t-1)(t+1))

Fit: r = 2, 3, 4.  Prediction verified at r = 5 (independent data).

This script:
  (1) prints the closed form
  (2) checks it against Rick's Day 198 data at r = 2, 3, 4, 5
  (3) prints the predicted value at r = 6, 7 (for eventual sober re-derivation).
"""
import sympy as sp

q, t = sp.symbols('q t')


def A_coef():
    return -(q**2 - 1) * (q - t - 1) / (t**2 - 1)

def B_coef():
    return t * (q**2 - 1) * (q - t) / (t - 1)

def C_coef():
    return -q * t**3 * (q**2 - 1) / (t**2 - 1)


def tau_closed(r_val):
    """The Day 200 closed form."""
    return (A_coef() + B_coef() * t**r_val + C_coef() * t**(2*r_val)) / q**3


def tau_data(r_val):
    """Day 198 SymPy-computed values."""
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
    raise ValueError(f"no data for r={r_val}")


def rewrite_closed_form():
    """
    Cast in a slightly cleaner form.  Note that
      A(q,t) = -(q^2-1)(q-t-1) / (t^2-1)
    We can pull out -(q^2-1)/(t^2-1) as a global factor.  Then

      tau_r * q^3 / [-(q^2-1)/(t^2-1)]
        = (q - t - 1) + (- t (q-t))(t+1) * t^r + q t^3 * t^(2r)
        = (q - t - 1) - t (q-t)(t+1) t^r + q t^(2r+3)

    Actually let me push the (t+1) into the B term properly:
      B = t(q^2-1)(q-t)/(t-1) = t(q^2-1)(q-t)(t+1) / (t^2-1).

    So writing K := -(q^2-1)/(t^2-1),
      tau_r * q^3 = K * [ (q-t-1) - (t+1) t (q-t) t^r + q t^3 * t^(2r) ]
                  = K * [ (q-t-1) - (t+1) (q-t) t^(r+1) + q t^(2r+3) ]

    So the closed form is:

      tau_r(q,t) = - (q^2-1) / [q^3 (t^2-1)]
                   * [ (q - t - 1)  -  (t+1)(q - t) t^(r+1)  +  q t^(2r+3) ]
    """
    K = -(q**2 - 1) / (t**2 - 1)
    for r_val in [2, 3, 4, 5]:
        pretty = K * ((q - t - 1) - (t + 1)*(q - t)*t**(r_val + 1) + q * t**(2*r_val + 3)) / q**3
        diff = sp.simplify(sp.expand(pretty - tau_data(r_val)))
        print(f"  r = {r_val}: pretty form matches data? {diff == 0}")
    return K


def qt1_specialization():
    """Sanity check at q = 1 (DS says tau_r|_{q=1} = 0)."""
    print()
    print("Specialization q = 1:  tau_r|_{q=1} =?= 0.")
    for r_val in [2, 3, 4, 5, 6, 7, 8]:
        val = sp.simplify(tau_closed(r_val).subs(q, 1))
        print(f"  r = {r_val}: tau_r(1,t) = {val}  (must be 0 by DS)")


def q0_expansion():
    """Also check q -> 0 behavior."""
    print()
    print("Behavior at q = 0:  tau_r(0, t) =")
    for r_val in [2, 3, 4, 5, 6, 7]:
        val = sp.series(tau_closed(r_val), q, 0, 1).removeO()
        # actually tau_r has 1/q^3 pole ==> Laurent
        # print numerator behavior for the leading term
        print(f"  r = {r_val}: leading q-behavior = {sp.factor(tau_closed(r_val) * q**3).subs(q, 0)}")


def report_all():
    print("=" * 72)
    print("Day 200 closed form for tau_r(q,t).")
    print("=" * 72)
    print()
    print("Coefficients (independent of r):")
    print(f"  A = {sp.factor(A_coef())}")
    print(f"  B = {sp.factor(B_coef())}")
    print(f"  C = {sp.factor(C_coef())}")
    print()
    print("Closed form:  tau_r(q,t) = (A + B t^r + C t^(2r)) / q^3.")
    print()
    print("Alternative pretty form (pull out K = -(q^2-1)/(t^2-1)):")
    print("  tau_r(q,t) = K/q^3 * [ (q-t-1) - (t+1)(q-t) t^(r+1) + q t^(2r+3) ]")
    print()
    print("Verification against Day 198 SymPy data:")
    for r_val in [2, 3, 4, 5]:
        pred = tau_closed(r_val)
        actual = tau_data(r_val)
        diff = sp.simplify(sp.expand(pred - actual))
        print(f"  r = {r_val}: closed - data = {diff}   {'OK' if diff == 0 else 'FAIL'}")

    rewrite_closed_form()
    qt1_specialization()

    print()
    print("Predictions (for eventual Day 201+ SymPy verification):")
    for r_val in [6, 7]:
        val = sp.factor(sp.simplify(tau_closed(r_val)))
        print(f"  r = {r_val}: tau_r = {val}")


if __name__ == "__main__":
    report_all()

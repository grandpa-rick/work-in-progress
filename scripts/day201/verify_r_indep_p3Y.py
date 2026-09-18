"""
Day 201: verify r-independent coefficients for p_3(Y) . e_r.

Conjectured r-indep formulas (based on r=2, r=3 data):
  c_{(r, 1, 1, 1)} = 1/q^6
  c_{(r, 2, 1)}    = -(q^2 t - q^2 + q t - q + t + 2)/q^6
  c_{(r+1, 1, 1)}  = (q^3 - 1)/q^6

More general conjecture for p_k(Y):
  c_{(r, 1^k)}         = 1 / q^{k(k+1)/2}
  c_{(r, 2, 1^{k-2})}  = -( q [k-1]_q (t-1) + (t + k - 1) ) / q^{k(k+1)/2}
  c_{(r+1, 1^{k-1})}   = (q^k - 1) / q^{k(k+1)/2}

Sanity check: at k=2, do these reduce to Rick's Lemma 1?
  c_{(r,1,1)} = 1/q^3               vs Lemma 1's 1/q^3.  MATCH.
  c_{(r,2)} = -(q(t-1) + t + 1)/q^3 = -(qt-q+t+1)/q^3  vs Lemma 1's -(qt-q+t+1)/q^3.  MATCH.
  c_{(r+1,1)} = (q^2-1)/q^3          vs Lemma 1's (q^2-1)/q^3.  MATCH.

We now test at k=3 against r=2, r=3 (already verified above).
Awaiting r=4 to confirm.
"""
import sympy as sp

q, t = sp.symbols('q t')

def q_int(n):
    """[n]_q = 1 + q + q^2 + ... + q^{n-1}."""
    return sum(q**i for i in range(n))

def leading_corner(k):
    """c_{(r, 1^k)} = 1/q^{k(k+1)/2}."""
    return sp.Rational(1) / q**(k*(k+1)//2)

def bump_first(k):
    """c_{(r+1, 1^{k-1})} = (q^k - 1)/q^{k(k+1)/2}."""
    return (q**k - 1) / q**(k*(k+1)//2)

def add_two_to_first(k):
    """c_{(r, 2, 1^{k-2})} = -(q [k-1]_q (t-1) + (t + k - 1))/q^{k(k+1)/2}."""
    return -(q * q_int(k-1) * (t-1) + (t + k - 1)) / q**(k*(k+1)//2)

def report(k):
    print(f"--- Conjectured r-independent coefficients for p_{k}(Y) ---")
    print(f"  c_(r, 1^{k})       = {sp.factor(leading_corner(k))}")
    print(f"  c_(r, 2, 1^{k-2})  = {sp.factor(add_two_to_first(k))}")
    print(f"  c_(r+1, 1^{k-1})   = {sp.factor(bump_first(k))}")

def check_k2():
    """Check Rick's Lemma 1 coefficients."""
    lemma1 = {
        (0, 1, 1): sp.Rational(1) / q**3,  # (r, 1, 1) -- pivot
        (0, 2): -(q*t - q + t + 1) / q**3,  # (r, 2)
        (1, 1): (q**2 - 1) / q**3,  # (r+1, 1)
    }
    print()
    print("k = 2 check vs Rick's Lemma 1:")
    diff1 = sp.simplify(leading_corner(2) - lemma1[(0, 1, 1)])
    diff2 = sp.simplify(add_two_to_first(2) - lemma1[(0, 2)])
    diff3 = sp.simplify(bump_first(2) - lemma1[(1, 1)])
    print(f"  c_(r,1,1) diff  = {diff1}")
    print(f"  c_(r,2) diff    = {diff2}")
    print(f"  c_(r+1,1) diff  = {diff3}")

def check_k3():
    """Check p_3(Y) coefficients at r=2, r=3."""
    # From p3Y_er.log
    r2 = {
        (0, 1, 1, 1): sp.Rational(1) / q**6,   # (2, 1, 1, 1)
        (0, 2, 1): -(q**2*t - q**2 + q*t - q + t + 2) / q**6,  # (2, 2, 1)
        (1, 1, 1): (q - 1)*(q**2 + q + 1) / q**6,  # (3, 1, 1)
    }
    print()
    print("k = 3 check vs r=2 data:")
    diff1 = sp.simplify(leading_corner(3) - r2[(0, 1, 1, 1)])
    diff2 = sp.simplify(add_two_to_first(3) - r2[(0, 2, 1)])
    diff3 = sp.simplify(bump_first(3) - r2[(1, 1, 1)])
    print(f"  c_(2,1,1,1) diff  = {diff1}")
    print(f"  c_(2,2,1) diff    = {diff2}")
    print(f"  c_(3,1,1) diff    = {diff3}")

if __name__ == "__main__":
    report(2)
    report(3)
    report(4)
    check_k2()
    check_k3()

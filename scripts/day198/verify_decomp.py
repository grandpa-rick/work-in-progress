"""Day 198 verification: C = p_2(Y).e_2 + 2t*D.

C = e_1 * e_1 * e_2  (Hikita star, Day 196 data)
D = e_2 * e_2        (Day 191)
p_2(Y).e_2           (Day 198 direct sub-agent computation)

Identity to check: for each partition mu of 4,
    [e_mu] C  ==  [e_mu] p_2(Y).e_2  +  2t * [e_mu] D
"""
import sympy as sp

q, t = sp.symbols('q t')

# Rick's Day 196 data (C = e_1 * e_1 * e_2):
c_C = {
    (4,):        (q - 1)**2 * (t + 1) * (t**2 + 1) * (q*t**2 + q*t + q + 1) / q**3,
    (3, 1):      (q - 1) * (2*q*t**2 + 2*q*t + q + 1) / q**3,
    (2, 2):      (q - 1) * (t + 1) / q**3,
    (2, 1, 1):   sp.Rational(1) / q**3,
    (1, 1, 1, 1): sp.Integer(0),
}

# Rick's Day 191 (D = e_2 * e_2):
c_D = {
    (4,):    (q - 1) * (1 + t**2) * (q*(1 + t + t**2) - t) / q**2,
    (3, 1):  (q - 1) * (t + 1) / q**2,
    (2, 2):  sp.Rational(1) / q**2,
    (2, 1, 1):    sp.Integer(0),
    (1, 1, 1, 1): sp.Integer(0),
}

# Sub-agent Day 198 direct computation (p_2(Y).e_2):
c_p2Ye2 = {
    (4,):        -(q - 1)*(q + 1)*(t**2 + 1)*(q*t**3 - q + t + 1) / q**3,
    (3, 1):       (q - 1)*(q + 1) / q**3,
    (2, 2):      -(q*t - q + t + 1) / q**3,
    (2, 1, 1):    sp.Rational(1) / q**3,
    (1, 1, 1, 1): sp.Integer(0),
}

print("Verifying  C  ==  p_2(Y).e_2  +  2t * D  coefficient-by-coefficient:\n")
all_ok = True
for mu in [(4,), (3,1), (2,2), (2,1,1), (1,1,1,1)]:
    lhs = sp.expand(c_C[mu])
    rhs = sp.expand(c_p2Ye2[mu] + 2*t*c_D[mu])
    diff = sp.simplify(sp.together(lhs - rhs))
    ok = (diff == 0)
    all_ok = all_ok and ok
    print(f"  mu = {mu}: diff = {diff}  {'OK' if ok else 'FAIL'}")

print()
print("ALL MATCH" if all_ok else "MISMATCH")

# Thm 4(1) symbolic check: zeta_free(ad_c(F_i)^n F_j) = t_i^n t_j q^{phi_ii C(n,2)} x^n prod_k (1 - q_ij q_ii^k y/x), all finite rank-2 types incl G2 n=4
import sympy as sp
exec(open('shuffle_checks.py').read().split("types = ")[0])
p11, p22, p12, p21, t1, t2 = sp.symbols('p11 p22 p12 p21 t1 t2', integer=True)
Q = sp.symbols('Q', positive=True)
types = {"A1xA1": ([[2, 0], [0, 2]], [1, 1]), "A2": ([[2, -1], [-1, 2]], [1, 1]), "B2": ([[2, -2], [-1, 2]], [1, 2]), "G2": ([[2, -3], [-1, 2]], [1, 3])}
P = [[p11, p12], [p21, p22]]; T = [t1, t2]
for name, (A, d) in types.items():
    Cs = [[d[i]*A[i][j] for j in range(2)] for i in range(2)]
    for i, j in [(0, 1), (1, 0)]:
        n = 1 - A[i][j]; S = adc_power(i, j, n, Cs, Q)
        lhs = 0
        for w, c in S.items():
            z = 1
            for a in range(len(w)):
                z *= T[w[a]]
                for b in range(a+1, len(w)): z *= Q**P[w[a]][w[b]]
            lhs += c*z
        x, y = Q**P[i][j], Q**P[j][i]
        rhs = T[i]**n*T[j]*Q**(P[i][i]*sp.binomial(n, 2))*x**n*sp.prod([1 - Q**Cs[i][j]*Q**(Cs[i][i]*k)*y/x for k in range(n)])
        print(name, (i, j), "n=", n, "Thm4(1) formula:", sp.simplify(sp.powsimp(sp.expand(lhs - rhs), force=True)) == 0)

"""Section 3(b) substitution check (Rick -> Clio, 2026-10-11).

Three independent computations of the two-row x two-part Green polynomial
X^{(l1,l2)}_{(x,y)}(t), defined by p_mu = sum_lambda X^lambda_mu(t) P_lambda(x;t):

 GT  ground truth: expand p_x p_y in TWO variables in the HL basis, where P_lambda(x1,x2;t)
     is built from Macdonald's symmetrisation definition III (2.1):
       P_lambda = (1/v_lambda(t)) sum_{w in S_2} w( x^lambda (x1 - t x2)/(x1 - x2) ).
     Two variables kill exactly the P_lambda with l(lambda)>2, so the coefficients of the
     two-row P_lambda are the two-row Green polynomials.  Shares nothing with Thm 2.5.
 S3b the Section 3(b) substitution (UID 784 PDF, 2026-10-07-reply-clio-III7-two-row.tex, Sec. 3 item (b)),
     evaluated term by term from G_{1,k}, then X = b_lambda * Phi_2 / ((1-t^x)(1-t^y)).
 YR  Young's-rule formula X = pi_k + (t-1) sum_{j<k} t^{k-1-j} pi_j, k = lambda_2,
     pi_j = #{j-subsets of [n] fixed by a permutation of cycle type (x,y)}.
 CF  closed form (FPSAC ex:tworow / 2026-10-08 cross-check note), for y <= x only.
"""
import sys
from itertools import combinations
import sympy as sp

t, x1, x2, w = sp.symbols('t x1 x2 w')
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10


def two_row_partitions(n):
    return [(n - k, k) for k in range(1, n // 2 + 1)]  # lambda_2 >= 1


def P2(lam):
    """HL P_lambda(x1,x2;t) by symmetrisation; v_lambda = 1 if l1>l2, (1+t) if l1=l2."""
    l1, l2 = lam
    f = x1**l1 * x2**l2 * (x1 - t * x2) / (x1 - x2)
    g = f + f.subs({x1: x2, x2: x1}, simultaneous=True)
    v = 1 if l1 > l2 else (1 + t)
    return sp.expand(sp.cancel(g / v))


def ground_truth(n, xx, yy):
    """Coefficients of P_lambda (all lambda with <=2 parts, |lambda|=n) in p_x p_y (x1,x2)."""
    f = sp.expand((x1**xx + x2**xx) * (x1**yy + x2**yy))
    lams = [(n, 0)] + two_row_partitions(n)
    # P_lambda = m_lambda + lower in dominance; peel from the top (lambda_1 largest first)
    coeffs = {}
    rem = f
    for lam in lams:  # (n,0),(n-1,1),... decreasing dominance
        c = sp.Poly(rem, x1, x2).coeff_monomial(x1**lam[0] * x2**lam[1])
        c = sp.factor(c)
        coeffs[lam] = sp.expand(c)
        Pl = P2(lam) if lam[1] > 0 else None
        if lam[1] == 0:  # P_(n)(x1,x2) by the same symmetrisation
            ff = x1**n * (x1 - t * x2) / (x1 - x2)
            Pl = sp.expand(sp.cancel(ff + ff.subs({x1: x2, x2: x1}, simultaneous=True)))
        rem = sp.expand(rem - c * Pl)
    assert sp.simplify(rem) == 0, rem
    return coeffs


def G1_coeffs(lam):
    """G_1(w) = P_rho(1,w;t), rho = lam - (1,1), as a dict k -> G_{1,k} (Sec. 3(b) display)."""
    rho1, rho2 = lam[0] - 1, lam[1] - 1
    m = rho1 - rho2
    if m == 0:
        return {rho2: sp.Integer(1)}
    d = {rho2: sp.Integer(1), rho2 + m: sp.Integer(1)}
    for k in range(1, m):
        d[rho2 + k] = 1 - t
    return d


def S3b(lam, xx, yy):
    """Phi_2 by the Sec. 3(b) formula, then X = b_lambda Phi_2/((1-t^x)(1-t^y))."""
    G = G1_coeffs(lam)
    g = lambda k: G.get(k, 0)
    rho1, rho2 = lam[0] - 1, lam[1] - 1
    Prho_1t = sum(c * t**k for k, c in G.items())  # P_rho(1,t;t) = G_1(t)
    inner = g(yy - 1) / (1 - t)**2 + sum((t**-j - t**j) * g(yy - 1 - j) for j in range(1, yy)) / (1 - t**2)
    brace = inner / t - Prho_1t / (1 - t**2) * (1 + t**-yy)
    Phi2 = (1 - t**xx) * (1 - t**yy) * brace
    b = (1 - t)**2 if lam[0] > lam[1] else (1 - t) * (1 - t**2)
    return sp.expand(sp.cancel(b * Phi2 / ((1 - t**xx) * (1 - t**yy)))), sp.factor(Phi2)


def young(lam, xx, yy):
    n = xx + yy
    k = lam[1]
    cyc = [set(range(xx)), set(range(xx, n))]
    def pi(j):
        return sum(1 for S in combinations(range(n), j)
                   if all(c <= set(S) or not (c & set(S)) for c in cyc))
    return sp.expand(pi(k) + (t - 1) * sum(t**(k - 1 - j) * pi(j) for j in range(k)))


def closed(lam, xx, yy):
    l1, l2 = lam
    assert yy <= xx
    mxy = 2 if xx == yy else 1
    if yy < l2:
        return sp.expand((t - 1) * t**(l2 - 1 - yy) * (1 + t**yy))
    if yy == l2:
        return sp.expand(mxy - (1 - t) * t**(l2 - 1))
    return sp.expand((t - 1) * t**(l2 - 1))


if __name__ == '__main__':
    ok = bad = 0
    rows = []
    for n in range(2, NMAX + 1):
        for xx in range(1, n):
            yy = n - xx
            gt = ground_truth(n, xx, yy)
            for lam in two_row_partitions(n):
                X_gt = gt[lam]
                X_s3b, _ = S3b(lam, xx, yy)
                X_yr = young(lam, xx, yy)
                checks = [X_s3b == X_gt, X_yr == X_gt]
                if yy <= xx:
                    checks.append(closed(lam, xx, yy) == X_gt)
                if all(checks):
                    ok += 1
                else:
                    bad += 1
                    print('FAIL', lam, (xx, yy), X_gt, X_s3b, X_yr)
                if yy <= xx and n <= 8:
                    rows.append((n, lam, (xx, yy), X_gt))
    print(f'ALL ordered (x,y) incl. y>x, |lambda|<={NMAX}: OK {ok} FAIL {bad}')
    print(f'y<=x rows with |lambda|<=8: {len(rows)}')
    # negative control: drop the j-sum (the (t^-j - t^j) term) -> must fail somewhere
    nc_fail = 0
    for n in range(2, 9):
        for xx in range(1, n):
            yy = n - xx
            gt = ground_truth(n, xx, yy)
            for lam in two_row_partitions(n):
                G = G1_coeffs(lam); g = lambda k: G.get(k, 0)
                Pr = sum(c * t**k for k, c in G.items())
                brace = (g(yy - 1) / (1 - t)**2) / t - Pr / (1 - t**2) * (1 + t**-yy)
                b = (1 - t)**2 if lam[0] > lam[1] else (1 - t) * (1 - t**2)
                if sp.expand(sp.cancel(b * brace)) != gt[lam]:
                    nc_fail += 1
    print(f'negative control (j-sum dropped), |lambda|<=8: {nc_fail} mismatches (must be >0)')
    # worked example
    lam, (xx, yy) = (5, 2), (4, 3)
    print('Worked example lambda=(5,2), (x,y)=(4,3):')
    print('  G_1 coeffs', G1_coeffs(lam))
    Xs, Phi = S3b(lam, xx, yy)
    print('  Phi_2 =', Phi, ' X(S3b) =', Xs, ' X(GT) =', ground_truth(7, 4, 3)[lam],
          ' X(Young) =', young(lam, xx, yy), ' closed =', closed(lam, xx, yy))
    print('TABLE')
    for n, lam, mu, X in rows:
        print(n, lam, mu, sp.factor(X), '|', X)

"""
Day 202: Symbolic proof-of-cancellation for Lemma 1 via R7 (Newton).

Given (as `computed`, r >= 2):
  Thm 3.12:  e_1 * e_r = alpha_r e_{r+1} + beta e_r e_1,
    alpha_r = (1 - q^-1)[r+1]_t,   beta = q^-1.
  Day 191:   e_2 * e_r = A_r e_{r+2} + B_r e_{r+1,1} + C_r e_{r,2},
    A_r = (1 - q^-1)[r+2]_t/[2]_t * ([r+1]_t - t[r-1]_t/q),
    B_r = q^-1 (1 - q^-1) [r]_t,
    C_r = q^-2.

R7 identity (analytic, proved via Newton in Lambda(Y) + intertwiner):
  p_2(Y).e_r  =  e_1 * (e_1 * e_r)  -  2t (e_2 * e_r).

We EXPAND e_1*(e_1*e_r) analytically in terms of atomic pieces:
  e_1 * (e_1*e_r)  =  alpha_r [e_1 * e_{r+1}]  +  beta [e_1 * (e_r e_1)].
  e_1 * e_{r+1}    =  alpha_{r+1} e_{r+2}  +  beta e_{r+1} e_1.
  e_1 * (e_r e_1)  =  ???                          <-- depth-2 obstruction
                     (this is Rick's "e_1 * (e_1 e_r)" tautology-producer).

BUT we can compute this depth-2 piece FROM THE DATA:
directly at r=2,3,4 SymPy gives its e-basis form; call it Z_r.

Then p_2(Y).e_r  =  alpha_r alpha_{r+1} e_{r+2}
                    + alpha_r beta e_{r+1,1}       (from e_1 * e_{r+1})
                    + beta * Z_r                     (from e_1 * (e_r e_1))
                    - 2t A_r e_{r+2}
                    - 2t B_r e_{r+1,1}
                    - 2t C_r e_{r,2}.

At the DS-relevant positions:
  (r,1,1):   0 + 0 + beta * c_(r,1,1)[Z_r] - 0 = q^-1 c_(r,1,1)[Z_r]
             = q^-3  (r-INDEP from Lemma 1).
  So c_(r,1,1)[Z_r] = q^-2.

  (r,2):     0 + 0 + beta * c_(r,2)[Z_r] - 2t C_r
             = q^-1 c_(r,2)[Z_r] - 2t q^-2
             = -(qt-q+t+1)/q^3
  So c_(r,2)[Z_r] = q^-2 (2t + (-qt+q-t-1)/(q))  ... let me just compute.

  (r+1,1):   alpha_r beta + beta c_(r+1,1)[Z_r] - 2t B_r
             = (1-q^-1)[r+1]_t / q + c_(r+1,1)[Z_r]/q - 2t q^-1(1-q^-1)[r]_t
             = (q^2-1)/q^3  (r-INDEP).

So c_(r+1,1)[Z_r] = q^2 [(q^2-1)/q^3 - (1-q^-1)[r+1]_t/q + 2t(1-q^-1)[r]_t/q ]
                  = (q^2-1)/q + q(q-1)/q * (-[r+1]_t + 2t [r]_t) / ... hmm messy.

The r-indep property of c_(r+1,1)[p_2(Y).e_r] = (q^2-1)/q^3 forces
   c_(r+1,1)[Z_r]  =  (q^2-1)/q  +  (1 - q^-1)([r+1]_t - 2t[r]_t) * q
(after multiplying by q).

Newton identity:  [r+1]_t - 2t[r]_t = 1 + t^{r+1} - 2t (1 + t + ... + t^{r-1})
                = 1 - t + t^{r+1} - 2 t^r + 2 t^r - ...  Actually
                [r+1]_t = 1 + t + t^2 + ... + t^r,  [r]_t = 1 + t + ... + t^{r-1}.
                [r+1]_t - 2t[r]_t = 1 + t + t^2 + ... + t^r  - 2t - 2t^2 - ... - 2t^r
                = 1 - t - t^2 - ... - t^r  = 1 - t[r]_t  = 1 - t(1-t^r)/(1-t).
Not obviously r-indep.

The point: Z_r itself IS r-dep (Rick's Day 191 flagged this).  But the
Newton-weighted linear combination
    beta * Z_r  -  2t * (e_2 * e_r)
kills the r-dep at (r+1, 1).

We verify this cancellation symbolically here.
"""
import sympy as sp

q, t = sp.symbols('q t')


def qint(n):
    if n <= 0:
        return sp.Integer(0)
    return sum(t**i for i in range(n))


# =======================================================================
# Compute Z_r := e_1 * (e_r * e_1)   [note e_r e_1 is the ORDINARY product]
# at rank m = r+2, expand in e-basis.
# =======================================================================

import sys
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day198')
from p2Y_er import build_action, e_r_X, expand_symmetric_in_e_basis, partitions_of


def part_sorted(lam):
    return tuple(sorted(lam, reverse=True))


def compute_Z_r(m, r):
    """
    Z_r := e_1 * (e_r(X) * e_1(X)) = e_1(Y) . (e_r(X) e_1(X)),
    since e_1(Y).F = e_1 * F for symmetric F.
    Degree = r + 2.
    """
    X, Ti_apply, Ti_inv_apply, Pi_apply, Y_apply = build_action(m)
    erX = e_r_X(m, r)
    e1X = e_r_X(m, 1)
    F = sp.expand(erX * e1X)
    total = sp.Integer(0)
    for i in range(1, m+1):
        total = sp.expand(total + Y_apply(F, i))
    return expand_symmetric_in_e_basis(total, m, r + 2)


def analytic_check(r):
    """
    Analytic decomposition:
      p_2(Y).e_r  =  alpha_r alpha_{r+1} e_{r+2}
                    + alpha_r beta e_{r+1,1}
                    + beta * Z_r
                    - 2t A_r e_{r+2}
                    - 2t B_r e_{r+1,1}
                    - 2t C_r e_{r,2}
    where
      alpha_r = (1 - q^-1)[r+1]_t
      beta = q^-1
      A_r = (1 - q^-1)[r+2]_t/[2]_t * ([r+1]_t - t[r-1]_t/q)
      B_r = q^-1 (1 - q^-1)[r]_t
      C_r = q^-2

    We evaluate the sub-leading positions and confirm r-independence.
    """
    alpha_r = (1 - 1/q) * qint(r + 1)
    alpha_rp1 = (1 - 1/q) * qint(r + 2)
    beta = 1/q
    A_r = (1 - 1/q) * qint(r + 2) / qint(2) * (qint(r + 1) - t * qint(r - 1) / q)
    B_r = (1/q) * (1 - 1/q) * qint(r)
    C_r = 1/q**2

    m = max(4, r + 2)
    Z_r = compute_Z_r(m, r)

    # Compute the analytic p_2(Y).e_r via R7 decomposition
    # We only track coefficients at the DS positions.
    positions = {
        (r + 2,):         (r+2,),
        part_sorted((r+1, 1)): 'r+1,1',
        part_sorted((r, 2)):   'r,2',
        part_sorted((r, 1, 1)):'r,1,1',
    }

    # Contribution at each position:
    # e_{r+2}:  alpha_r alpha_{r+1}  +  beta * c[Z_r,(r+2)]  -  2t A_r
    # e_{r+1,1}: alpha_r beta + beta * c[Z_r,(r+1,1)] - 2t B_r
    # e_{r,2}:  beta * c[Z_r,(r,2)] - 2t C_r
    # e_{r,1,1}: beta * c[Z_r,(r,1,1)]

    def c(D, lam):
        return sp.simplify(D.get(lam, sp.Integer(0)))

    contribs = {}
    contribs[(r+2,)] = sp.simplify(alpha_r * alpha_rp1 + beta * c(Z_r, (r+2,)) - 2*t*A_r)
    contribs[part_sorted((r+1,1))] = sp.simplify(
        alpha_r * beta + beta * c(Z_r, part_sorted((r+1,1))) - 2*t*B_r)
    contribs[part_sorted((r,2))] = sp.simplify(
        beta * c(Z_r, part_sorted((r,2))) - 2*t*C_r)
    contribs[part_sorted((r,1,1))] = sp.simplify(
        beta * c(Z_r, part_sorted((r,1,1))))

    print(f"\n=== Analytic R7 check at r={r} ===")
    for lam, label in positions.items():
        val = sp.factor(contribs.get(lam, sp.Integer(0)))
        print(f"  c_{lam}[p_2(Y).e_r] (analytic) = {val}")

    # Print Z_r itself for reference
    print(f"\n  Z_r = e_1 * (e_r e_1) individual coefficients:")
    for lam in partitions_of(r + 2):
        val = sp.factor(sp.simplify(Z_r.get(lam, sp.Integer(0))))
        if val != 0:
            print(f"    c_{lam}[Z_r] = {val}")


# Now the actual test of r-independence: compare r=2, r=3, r=4 sub-leading positions
def r_indep_check():
    print("\n\n" + "=" * 76)
    print("Direct r-independence proof from R7 decomposition")
    print("=" * 76)

    # For r=2, 3, 4 compute the analytic sub-leading coefficients
    # via the R7 recipe, and check that (r,1,1), (r,2), (r+1,1) match
    # position-shifted between r's.
    results = {}
    for r in [2, 3, 4, 5]:
        m = max(4, r + 2)
        Z_r = compute_Z_r(m, r)
        alpha_r = (1 - 1/q) * qint(r + 1)
        beta = 1/q
        B_r = (1/q) * (1 - 1/q) * qint(r)
        C_r_coef = 1/q**2

        # (r, 1, 1):  = beta * c_(r,1,1)[Z_r]
        c_r11 = beta * sp.simplify(Z_r.get(part_sorted((r, 1, 1)), sp.Integer(0)))
        # (r, 2):    = beta * c_(r,2)[Z_r] - 2 t C_r
        c_r2 = beta * sp.simplify(Z_r.get(part_sorted((r, 2)), sp.Integer(0))) - 2*t*C_r_coef
        # (r+1, 1):  = alpha_r * beta + beta * c_(r+1,1)[Z_r] - 2 t B_r
        c_rp11 = alpha_r*beta + beta * sp.simplify(Z_r.get(part_sorted((r+1, 1)), sp.Integer(0))) - 2*t*B_r
        # (r+2):     = alpha_r * alpha_{r+1} + beta * c_(r+2)[Z_r] - 2 t A_r
        alpha_rp1 = (1 - 1/q) * qint(r + 2)
        A_r_val = (1 - 1/q) * qint(r + 2) / qint(2) * (qint(r + 1) - t * qint(r - 1) / q)
        c_rp2 = alpha_r*alpha_rp1 + beta * sp.simplify(Z_r.get((r+2,), sp.Integer(0))) - 2*t*A_r_val

        results[r] = {
            (r, 1, 1): sp.factor(sp.simplify(c_r11)),
            (r, 2):    sp.factor(sp.simplify(c_r2)),
            (r+1, 1):  sp.factor(sp.simplify(c_rp11)),
            (r+2,):    sp.factor(sp.simplify(c_rp2)),
        }
        print(f"\n  r={r}:")
        for pos, val in results[r].items():
            print(f"    c_{pos} = {val}")

    # r-indep test: compare r=2 vs r=3 vs r=4 vs r=5 at position-labels (r,1,1),(r,2),(r+1,1)
    print("\n  R-INDEP TEST: subtract successive r's for each sub-leading position.")
    for pos_name in ['(r,1,1)', '(r,2)', '(r+1,1)']:
        vals = []
        for r in [2, 3, 4, 5]:
            if pos_name == '(r,1,1)':
                p = (r, 1, 1)
            elif pos_name == '(r,2)':
                p = (r, 2)
            else:
                p = (r+1, 1)
            vals.append(results[r][p])
        print(f"    Position {pos_name}: values across r=2..5: {vals}")
        diffs = [sp.simplify(vals[i+1] - vals[i]) for i in range(len(vals)-1)]
        print(f"      differences: {diffs}   [should all be 0 for r-indep]")


if __name__ == "__main__":
    for r in [2, 3, 4]:
        analytic_check(r)
    r_indep_check()

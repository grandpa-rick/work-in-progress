"""
Day 204: Reduction of Sub-Lemma Z to four partial-symmetrizer identities (L1)-(L4).

Sub-Lemma Z states, for r >= 2:
    Z_r := e_1 * e_(r,1) = e_1(Y) . (e_r(X) * e_1(X))
has four nonzero e-basis coefficients:
    c_(r,1,1) = 1/q^2
    c_(r,2)   = (q-1)(t+1)/q^2
    c_(r+1,1) = (q-1)(q t [r]_t + 1)/q^2
    c_(r+2)   = (q-1)^2 [r+2]_t/q^2

STRATEGY:
    Step A. e_1(Y) . F = sigma_m . pi . F, where sigma_m = sum_{k=0..m-1} T_k T_(k-1) ... T_1
            (identity on symmetric polynomials of degree <= m). VERIFIED HERE.

    Step B. pi(F G) = pi(F) pi(G) / X_1. Applied to F=e_r, G=e_1:
            pi(e_r * e_1) = X_1 . fh + q^-1 X_1^2 . (f + gh) + q^-2 X_1^3 . g,
            where f = e_r(tail), g = e_(r-1)(tail), h = e_1(tail).

    Step C. Apply sigma_m to each piece separately (piece is X_1^a . G_tail).

    Step D. Combine and match Sub-Lemma Z's four coefficients.

RESIDUAL GAP: prove (L1)-(L4) symbolic in r (all four are q-free identities in the
AHA level-1 rep).
"""
import sympy as sp
from itertools import combinations
import sys

sys.path.insert(0, '/home/agent/projects/proofs/scripts/day198')
from p2Y_er import (
    build_action, e_r_X, expand_symmetric_in_e_basis, partitions_of,
)

q, t = sp.symbols('q t')


def qint(n):
    if n <= 0:
        return sp.Integer(0)
    return sum(t**i for i in range(n))


def part_sorted(lam):
    return tuple(sorted(lam, reverse=True))


def e_tail(m, r):
    """e_r(X_2, ..., X_m)."""
    if r == 0:
        return sp.Integer(1)
    Xs = sp.symbols(f'X1:{m+1}')
    Xtail = Xs[1:]
    if r > len(Xtail):
        return sp.Integer(0)
    return sum(sp.prod(c) for c in combinations(Xtail, r))


def sigma_m_apply(F, m):
    """Apply sigma_m = sum_{k=0..m-1} T_k T_(k-1) ... T_1 to F.
    (k=0 term is identity.)"""
    _, Ti_apply, _, _, _ = build_action(m)
    total = F
    for k in range(1, m):
        curr = F
        for j in range(1, k + 1):
            curr = Ti_apply(curr, j)
        total = sp.expand(total + curr)
    return total


def sigma_m_pi(F, m):
    """sigma_m . pi . F."""
    _, _, _, Pi_apply, _ = build_action(m)
    return sigma_m_apply(Pi_apply(F), m)


# ================================================================
# Step A: verify e_1(Y) . F = sigma_m . pi . F for symmetric F.
# ================================================================
def verify_sigma_m_pi_equals_e1_Y():
    print("=" * 76)
    print("Step A: Verify e_1(Y) . F = sigma_m . pi . F for symmetric F")
    print("=" * 76)
    all_ok = True
    for m in [3, 4, 5]:
        _, _, _, _, Y_apply = build_action(m)
        tests = [
            ('e_1', e_r_X(m, 1)),
            ('e_2', e_r_X(m, 2)),
            ('e_1^2', sp.expand(e_r_X(m, 1)**2)),
            ('e_1 * e_2', sp.expand(e_r_X(m, 1) * e_r_X(m, 2))),
            ('e_1 * e_3', sp.expand(e_r_X(m, 1) * e_r_X(m, 3))),
        ]
        for name, F in tests:
            if F == 0:
                continue
            lhs = sigma_m_pi(F, m)
            rhs = sum(Y_apply(F, i) for i in range(1, m + 1))
            d = sp.simplify(sp.expand(lhs - rhs))
            status = "OK" if d == 0 else "FAIL"
            if d != 0:
                all_ok = False
            print(f"  m={m}, F={name}: diff = {d}  [{status}]")
    print(f"\nStep A overall: {'PASSED' if all_ok else 'FAILED'}")
    return all_ok


# ================================================================
# Step B: verify pi decomposition of e_r * e_1.
# ================================================================
def verify_pi_split():
    print("\n" + "=" * 76)
    print("Step B: Verify pi(e_r * e_1) = X_1 f h + q^-1 X_1^2 (f + g h) + q^-2 X_1^3 g")
    print("=" * 76)
    all_ok = True
    for r in [1, 2, 3, 4]:
        m = r + 2
        _, _, _, Pi_apply, _ = build_action(m)
        Xs = sp.symbols(f'X1:{m+1}')
        f = e_tail(m, r)
        g = e_tail(m, r - 1)
        h = e_tail(m, 1)
        lhs = sp.expand(Pi_apply(sp.expand(e_r_X(m, r) * e_r_X(m, 1))))
        rhs = sp.expand(Xs[0] * f * h
                        + Xs[0]**2 / q * (f + g * h)
                        + Xs[0]**3 / q**2 * g)
        d = sp.simplify(lhs - rhs)
        status = "OK" if d == 0 else "FAIL"
        if d != 0:
            all_ok = False
        print(f"  r={r}, m={m}: pi(e_r * e_1) - decomposition = {d}  [{status}]")
    print(f"\nStep B overall: {'PASSED' if all_ok else 'FAILED'}")
    return all_ok


# ================================================================
# Step C: compute each of the four pieces via sigma_m and expand.
# ================================================================
def compute_piece(m, r, piece_id):
    """
    piece_id = 'L1': sigma_m . X_1 * e_r(tail) e_1(tail)
             = 'L2': sigma_m . X_1^2 * e_r(tail)
             = 'L3': sigma_m . X_1^2 * e_(r-1)(tail) e_1(tail)
             = 'L4': sigma_m . X_1^3 * e_(r-1)(tail)
    Returns dict of e-basis coefficients (partitions of r+2).
    """
    Xs = sp.symbols(f'X1:{m+1}')
    if piece_id == 'L1':
        F = sp.expand(Xs[0] * e_tail(m, r) * e_tail(m, 1))
    elif piece_id == 'L2':
        F = sp.expand(Xs[0]**2 * e_tail(m, r))
    elif piece_id == 'L3':
        F = sp.expand(Xs[0]**2 * e_tail(m, r - 1) * e_tail(m, 1))
    elif piece_id == 'L4':
        F = sp.expand(Xs[0]**3 * e_tail(m, r - 1))
    else:
        raise ValueError(piece_id)
    val = sigma_m_apply(F, m)
    return expand_symmetric_in_e_basis(val, m, r + 2)


def claimed_piece(r, piece_id):
    """The boxed formula for each of (L1)-(L4)."""
    if piece_id == 'L1':
        return {
            (r + 2,): qint(r + 2),
            part_sorted((r + 1, 1)): t * qint(r),
        }
    if piece_id == 'L2':
        return {
            (r + 2,): -qint(r + 2),
            part_sorted((r + 1, 1)): sp.Integer(1),
        }
    if piece_id == 'L3':
        return {
            (r + 2,): -qint(r + 2),
            part_sorted((r + 1, 1)): -t * qint(r),
            part_sorted((r, 2)): qint(2),
        }
    if piece_id == 'L4':
        return {
            (r + 2,): qint(r + 2),
            part_sorted((r + 1, 1)): -sp.Integer(1),
            part_sorted((r, 2)): -qint(2),
            part_sorted((r, 1, 1)): sp.Integer(1),
        }
    raise ValueError(piece_id)


def verify_pieces():
    print("\n" + "=" * 76)
    print("Step C: Verify (L1)-(L4) at r = 2, 3, 4, 5 (r=1 has partition collision)")
    print("=" * 76)
    all_ok = True
    # Note: at r = 1, partitions (r,2) = (r+1,1) = (2,1) collide, so the claimed_piece
    # entries have to be *summed* at the collision. Sub-Lemma Z is stated for r >= 2.
    for r in [2, 3, 4, 5]:
        m = r + 2
        print(f"\n--- r = {r}, m = {m} ---")
        for pid in ['L1', 'L2', 'L3', 'L4']:
            actual = compute_piece(m, r, pid)
            claim = claimed_piece(r, pid)
            piece_ok = True
            for lam in partitions_of(r + 2):
                a_val = sp.simplify(actual.get(lam, sp.Integer(0)))
                c_val = sp.simplify(claim.get(lam, sp.Integer(0)))
                diff = sp.simplify(sp.together(a_val - c_val))
                if diff != 0:
                    piece_ok = False
                    all_ok = False
                    print(f"  {pid}, {lam}: actual = {sp.factor(a_val)}, claimed = {sp.factor(c_val)}  [FAIL]")
            status = "OK" if piece_ok else "FAIL"
            print(f"  {pid}: {status}")
    print(f"\nStep C overall: {'PASSED' if all_ok else 'FAILED'}")
    return all_ok


# ================================================================
# Step D: assemble (L1)-(L4) with weights (1, 1/q, 1/q, 1/q^2) and
# match against Sub-Lemma Z.
# ================================================================
def claimed_Z(r):
    return {
        part_sorted((r, 1, 1)): sp.Rational(1) / q**2,
        part_sorted((r, 2)): (q - 1) * (t + 1) / q**2,
        part_sorted((r + 1, 1)): (q - 1) * (q * t * qint(r) + 1) / q**2,
        (r + 2,): (q - 1)**2 * qint(r + 2) / q**2,
    }


def verify_assembly():
    print("\n" + "=" * 76)
    print("Step D: assemble (L1)-(L4) -> match Sub-Lemma Z at r = 2, 3, 4, 5")
    print("=" * 76)
    all_ok = True
    weights = {'L1': sp.Integer(1), 'L2': sp.Rational(1) / q,
               'L3': sp.Rational(1) / q, 'L4': sp.Rational(1) / q**2}
    for r in [2, 3, 4, 5]:
        m = r + 2
        assembled = {lam: sp.Integer(0) for lam in partitions_of(r + 2)}
        for pid in ['L1', 'L2', 'L3', 'L4']:
            claim = claimed_piece(r, pid)
            for lam, c in claim.items():
                assembled[lam] = sp.expand(assembled[lam] + weights[pid] * c)
        target = claimed_Z(r)
        print(f"\n--- r = {r} ---")
        for lam in partitions_of(r + 2):
            a_val = sp.factor(sp.simplify(assembled[lam]))
            t_val = sp.factor(sp.simplify(target.get(lam, sp.Integer(0))))
            diff = sp.simplify(sp.together(a_val - t_val))
            status = "OK" if diff == 0 else "FAIL"
            if diff != 0:
                all_ok = False
            marker = "" if lam in target else "  (should be 0)"
            print(f"  {lam}: assembled = {a_val}, target = {t_val}{marker}  [{status}]")
    print(f"\nStep D overall: {'PASSED' if all_ok else 'FAILED'}")
    return all_ok


def main():
    print()
    print("Day 204 attempt: Sub-Lemma Z via sigma_m . pi decomposition")
    print()
    ok_a = verify_sigma_m_pi_equals_e1_Y()
    ok_b = verify_pi_split()
    ok_c = verify_pieces()
    ok_d = verify_assembly()
    print("\n" + "=" * 76)
    print(f"OVERALL: A={ok_a}, B={ok_b}, C={ok_c}, D={ok_d}")
    if ok_a and ok_b and ok_c and ok_d:
        print("All four steps verified computationally. Analytic gap: prove (L1)-(L4)")
        print("symbolic in r (q-free identities in the AHA level-1 rep).")
    print("=" * 76)


if __name__ == "__main__":
    main()

"""
Day 207: cross-check of DS(r,1,1) direct route.

Computes C_r = e_1(Y) . (e_1(Y) . e_r)  [= e_1 * e_1 * e_r by R0]
directly in the AHA level-1 polynomial rep (day205/k3_fast_pipeline.AHA,
conventions = day198/p2Y_er.py::build_action), expands in the e-basis,
and compares with the closed form (s = 1/q, [n] = [n]_t):

  c_(r+2)    = (1-s)^2 [r+2] ([r+1] + s)
  c_(r+1,1)  = s (1-s) ([r+1] + t[r] + s)
  c_(r,2)    = s^2 (1-s) [2]
  c_(r,1,1)  = s^3
  all other mu: 0

Also checks the intermediate e_1(Y).e_{r+1} (Thm 3.12) and e_1(Y).(e_1 e_r) (Sub-Lemma Z).
Usage: python3 ds_r11_direct_check.py
"""
import sys, time
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
import sympy as sp
from k3_fast_pipeline import AHA, e_expansion, to_sympy

q, t = sp.symbols('q t')
s = 1 / q


def qi(n):
    return sum(t**i for i in range(n)) if n > 0 else 0


def e1Y(A, F):
    tot = 0 * A.s
    for i in range(1, A.m + 1):
        tot = tot + A.Y(F, i)
    return tot


def check(expanded, claim, label):
    keys = set(expanded) | set(claim)
    ok = True
    for mu in keys:
        d = sp.simplify(expanded.get(mu, 0) - claim.get(mu, 0))
        if d != 0:
            ok = False
            print(f'   MISMATCH {label} {mu}: diff={sp.factor(d)}')
    return ok


def main():
    allok = True
    for r in range(2, 6):
        for m in (r + 2, r + 3):
            if r == 5 and m == 8:
                pass
            t0 = time.time()
            A = AHA(m)
            er, er1, e1 = A.e(r), A.e(r + 1), A.e(1)
            # Thm 3.12 at r+1
            thm = {lam: to_sympy(A, c) for lam, c in e_expansion(A, e1Y(A, er1), r + 2).items()}
            thm_claim = {(r + 2,): (1 - s) * qi(r + 2), (r + 1, 1): s}
            ok1 = check(thm, thm_claim, 'Thm3.12')
            # Sub-Lemma Z
            Z = {lam: to_sympy(A, c) for lam, c in e_expansion(A, e1Y(A, e1 * er), r + 2).items()}
            Z_claim = {(r + 2,): (1 - s)**2 * qi(r + 2), (r + 1, 1): (1 - s) * (t * qi(r) + s),
                       (r, 2): s * (1 - s) * qi(2), (r, 1, 1): s**2}
            ok2 = check(Z, Z_claim, 'SubLemmaZ')
            # C_r directly
            C = {lam: to_sympy(A, c) for lam, c in e_expansion(A, e1Y(A, e1Y(A, er)), r + 2).items()}
            C_claim = {(r + 2,): (1 - s)**2 * qi(r + 2) * (qi(r + 1) + s),
                       (r + 1, 1): s * (1 - s) * (qi(r + 1) + t * qi(r) + s),
                       (r, 2): s**2 * (1 - s) * qi(2), (r, 1, 1): s**3}
            ok3 = check(C, C_claim, 'C_r')
            supp = sorted(mu for mu, c in C.items() if sp.simplify(c) != 0)
            ok = ok1 and ok2 and ok3
            allok &= ok
            print(f'r={r} m={m}: Thm3.12 {ok1}  Z {ok2}  C_r {ok3}  support={supp}  '
                  f'lead={sp.simplify(C.get((r,1,1),0))}  ({time.time()-t0:.1f}s)', flush=True)
    print('ALL OK' if allok else 'FAILURES')


if __name__ == '__main__':
    main()

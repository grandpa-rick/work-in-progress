# Day 243 (2026-10-17): two controls, two blocks (Clio review 2026-10-16 §3).
# Arm A: inverted braiding in the bialgebra middle-crossing -> must break bialg, cannot break assoc.
# Arm B: non-bilinear phi (fails 2-cocycle) -> must break assoc, cannot break bialg.
# Same domain as lattice_bialg.py.
from lattice import *
from lattice_bialg import comps2, Delta, tens
class QQnb(QQ):
    def ph(s, a, b): return s.q ** (a[0]*b[0]*(a[1]+b[1]+1) + a[1]*b[1]*b[1])  # not bilinear
def run(Q, csgn):
    bad = tot = ba = bt = 0
    degs = [(1, 0), (0, 1), (1, 1), (2, 0), (2, 1)]
    for d1 in degs:
        for d2 in degs:
            for a in comps2(d1):
                for b in comps2(d2):
                    L = {k: v for k, v in Delta(Q.mul({a: Fr(1)}, {b: Fr(1)})).items() if v != 0}
                    R = tens(Q, Delta({a: Fr(1)}), Delta({b: Fr(1)}), csgn)
                    tot += 1; bad += L != R
                    for d3 in degs[:3]:
                        for c in comps2(d3):
                            A_, B_, C_ = {a: Fr(1)}, {b: Fr(1)}, {c: Fr(1)}
                            bt += 1; ba += Q.mul(Q.mul(A_, B_), C_) != Q.mul(A_, Q.mul(B_, C_))
    return bad, tot, ba, bt
q = Fr(3, 7); C = [[2, -1], [-1, 2]]
for name, Q, s in [("bilinear phi, true braiding", QQ(C, [[1, -1], [0, 1]], q), 1),
                   ("ARM A inverted braiding", QQ(C, [[1, -1], [0, 1]], q), -1),
                   ("ARM B non-bilinear phi", QQnb(C, None, q), 1)]:
    bad, tot, ba, bt = run(Q, s)
    print(f"{name:30s} bialg FAIL {bad}/{tot}   assoc FAIL {ba}/{bt}")

# Day 242 PROVE: A_3 separator for Thm 4(4) + sorted-word character test.
# 'f = free braided algebra on F_0,F_1,F_2; Theta = quantum shuffle map (simple-letter part of Psi, t=1);
# f_nu := 'f_nu / ker Theta. Psi_phi computed on 'f via Lemma A (lattice.QQ.Psi).
import sys
from lattice import *
from itertools import permutations
q = Fr(int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else Fr(3, 7)
C = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
I = 3
t = [Fr(2), Fr(-5, 3), Fr(7, 4)]
def words_of(nu):
    base = [i for i in range(I) for _ in range(nu[i])]
    return sorted(set(permutations(base)))
def phi_order(order):  # phi(a_k,a_l) = (a_k,a_l) if k after l in order, else 0 (diag 0)
    pos = {v: n for n, v in enumerate(order)}
    return [[C[k][l] if pos[k] > pos[l] else 0 for l in range(I)] for k in range(I)]
def theta(Q, w):
    P = Q.Psi(w, [Fr(1)]*I)
    return {k: v for k, v in P.items() if all(sum(x) == 1 for x in k)}
def mat(cols):
    keys = sorted({k for c in cols for k in c})
    return sp.Matrix([[sp.Rational(c.get(k, 0).numerator, c.get(k, 0).denominator) if k in c else 0 for c in cols] for k in keys]) if keys else sp.zeros(1, len(cols))
nus = [(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1), (1, 2, 2)]
Bs = {"O1 order 0<1<2 (1->2->3)": phi_order([0, 1, 2]),
      "O2 order 0<2<1 (1->2<-3)": phi_order([0, 2, 1]),
      "O1+random symmetric": [[phi_order([0, 1, 2])[i][j] + S for j, S in enumerate(row)] for i, row in enumerate([[1, 2, -1], [2, 0, 3], [-1, 3, -2]])],
      "WRONG antisym 0 on edge 1-2": [[0, 0, 0], [-1, 0, 0], [0, 0, 0]],
      "phi=0": [[0]*3 for _ in range(3)]}
Q0 = QQ(C, [[0]*3 for _ in range(3)], q)
for nu in nus:
    W = words_of(nu)
    T = mat([theta(Q0, w) for w in W])
    K = T.nullspace()  # ker Theta in 'f_nu (coefficient vectors over words)
    dimf = len(W) - len(K)
    print(f"nu={nu}: #words={len(W)} dim f_nu={dimf} dim kerTheta={len(K)}")
    for lab, B in Bs.items():
        Q = QQ(C, B, q)
        cols = [Q.Psi(w, t) for w in W]
        M = mat(cols)
        r = M.rank()
        desc = all((M*k).is_zero_matrix for k in K)
        print(f"   {lab:32s} rank Psi on 'f_nu={r}  descends to f (Psi(kerTheta)=0)? {desc}  injective on f_nu? {desc and r == dimf}")

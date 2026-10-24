# Wake 255 / Q34 follow-up: T is multiplicative over components (exp-type); MN for p_r is connected (log-type).
# So also compare G := log T = (1/2) log(I+C) with c_r * [Q_kappa] p_r Q_rho, r = |kappa/rho|.
import sys, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
sys.argv = [sys.argv[0], str(N)]
sys.path.insert(0, os.path.join(HERE, '../prove254'))
cwd = os.getcwd(); os.chdir(os.path.join(HERE, '../prove254')); from half import KS, SUP; os.chdir(cwd)
from qmn import MN, boxes
idx = {k: i for i, k in enumerate(KS)}; n = len(KS)
C = [[Fr(0)]*n for _ in range(n)]
for k, d in SUP.items():
    for r, v in d.items(): C[idx[k]][idx[r]] = Fr(v)
def mul(A, B): return [[sum(A[i][l]*B[l][j] for l in range(n) if A[i][l]) for j in range(n)] for i in range(n)]
G = [[Fr(0)]*n for _ in range(n)]; P = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
for e in range(1, sum(KS[-1]) + 2):
    P = mul(P, C); co = Fr((-1)**(e+1), 2*e)
    G = [[G[i][j] + co*P[i][j] for j in range(n)] for i in range(n)]
from collections import defaultdict
by = defaultdict(list)
for k in KS:
    for r in KS:
        if boxes(r) < boxes(k):
            s = sum(k) - sum(r); g = G[idx[k]][idx[r]]
            m = MN(r, s).get(k, Fr(0)) if s % 2 else None
            by[s].append((k, r, g, m))
for s in sorted(by):
    L = by[s]; nzG = [x for x in L if x[2]]
    if s % 2:
        rat = {str(x[2]/x[3]) if x[3] else ('G/0' if x[2] else '0/0') for x in L}
        print(f'r={s}: #G!=0 {len(nzG)}, #MN!=0 {sum(1 for x in L if x[3])}, ratios G/MN {sorted(rat)[:10]}')
    else:
        print(f'r={s} (even, no p_r in Gamma): #G!=0 {len(nzG)} of {len(L)}; e.g. {[(x[0],x[1],str(x[2])) for x in nzG[:3]]}')

# Wake 255 / Q34: compare T=(I+C)^{1/2} entries with shifted MN coefficients [Q_kappa] p_r Q_rho (r odd),
# and test the more general hunch "T = multiplication by one element g of (completed) Gamma in the Q-basis".
import sys, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
sys.argv = [sys.argv[0], str(N)]
sys.path.insert(0, os.path.join(HERE, '../prove254'))
cwd = os.getcwd(); os.chdir(os.path.join(HERE, '../prove254'))
from half import KS, SUP          # strict kappa with parts <= N, |c| = V(theta) (Q30)
os.chdir(cwd)
from qmn import MN, Q, expand_Q, boxes, mul_p
idx = {k: i for i, k in enumerate(KS)}; n = len(KS)
C = [[Fr(0)]*n for _ in range(n)]
for k, d in SUP.items():
    for r, v in d.items(): C[idx[k]][idx[r]] = Fr(v)
def mul(A, B): return [[sum(A[i][l]*B[l][j] for l in range(n) if A[i][l]) for j in range(n)] for i in range(n)]
T = [[Fr(int(i == j)) for j in range(n)] for i in range(n)]; P = [row[:] for row in T]; b = Fr(1)
for e in range(1, sum(KS[-1]) + 2):   # C strictly lowers |kappa| => nilpotent
    P = mul(P, C); b = b * (Fr(1, 2) - (e - 1)) / e
    T = [[T[i][j] + b*P[i][j] for j in range(n)] for i in range(n)]
def Tv(k, r): return T[idx[k]][idx[r]]
pairs = [(k, r) for k in KS for r in KS if boxes(r) < boxes(k)]
# (1) direct MN comparison on odd sizes
rows = []; agree = dis = 0; firstdis = None
for k, r in sorted(pairs, key=lambda p: (sum(p[0]) - sum(p[1]), sum(p[0]))):
    s = sum(k) - sum(r)
    if s % 2 == 0: continue
    m = MN(r, s).get(k, Fr(0)); t = Tv(k, r)
    if m == t: agree += 1
    else:
        dis += 1
        if firstdis is None: firstdis = (k, r, s, t, m)
    rows.append((k, r, s, t, m))
print(f'N={N}: odd-size pairs {agree+dis}: T == [Q_k] p_r Q_rho on {agree}, differ on {dis}')
print('first discrepancy (kappa, rho, r, T, MN):', firstdis)
# ratio test per r: is T_r = c_r * MN_r ?
for s in sorted({x[2] for x in rows}):
    rat = {(x[3] / x[4]) if x[4] else ('T/0' if x[3] else '0/0') for x in rows if x[2] == s}
    supT = sum(1 for x in rows if x[2] == s and x[3]); supM = sum(1 for x in rows if x[2] == s and x[4])
    both = sum(1 for x in rows if x[2] == s and x[3] and x[4])
    print(f'  r={s}: #T!=0 {supT}, #MN!=0 {supM}, both {both}, ratios T/MN {sorted(map(str,rat))[:8]}')
print('sample (kappa, rho, r, T, MN):')
for x in rows[:40:3]: print('  ', x[0], x[1], x[2], str(x[3]), str(x[4]))
# (2) multiplication hypothesis: g := sum_kappa T(kappa, 0) Q_kappa; is T(kappa,rho) = [Q_kappa] g Q_rho ?
g = {k: Tv(k, ()) for k in KS if Tv(k, ())}
print('g = T(.,0) in Q-basis:', {k: str(v) for k, v in g.items()})
def prod_Q(a, b, nv):  # [Q] coefficients of Q_a Q_b
    qa, qb = Q(a, nv), Q(b, nv); out = {}
    for e1, c1 in qa.items():
        for e2, c2 in qb.items():
            e = tuple(x + y for x, y in zip(e1, e2)); out[e] = out.get(e, 0) + c1*c2
    return expand_Q(out, nv)
ok = bad = 0; fb = None; maxdeg = sum(KS[-1])
for k, r in pairs:
    deg = sum(k)
    pred = sum((v * prod_Q(a, r, deg).get(k, 0) for a, v in g.items() if sum(a) + sum(r) == deg), Fr(0))
    if pred == Tv(k, r): ok += 1
    else:
        bad += 1
        if fb is None or (sum(k)-sum(r), sum(k)) < (sum(fb[0])-sum(fb[1]), sum(fb[0])): fb = (k, r, Tv(k, r), pred)
print(f'multiplication-by-g test: {ok} ok, {bad} fail; smallest failure (kappa, rho, T, [Q_k] g Q_rho): {fb}')

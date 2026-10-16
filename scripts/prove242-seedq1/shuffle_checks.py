# Day 242: (1) Serre elements lie in ker Theta (quantum shuffle), symbolic q, all finite rank-2 types.
#          (2) sorted-word functional lambda_o is a phi_o-character on Sh_q(V), and zeta_free = lambda_o o Theta o D_t on 'f.
import sympy as sp, random
from itertools import product as iprod
from collections import defaultdict
q = sp.symbols('q')
def shmul(u, v, C, qq):  # quantum shuffle of words u,v: v-letter placed before u-letter gives qq^{(v,u)}
    memo = {}
    def rec(i, j):
        if (i, j) in memo: return memo[(i, j)]
        out = defaultdict(lambda: 0)
        if i == len(u) and j == len(v): out[()] = 1
        if i < len(u):
            for w, c in rec(i+1, j).items(): out[(u[i],)+w] += c
        if j < len(v):
            wt = qq**sum(C[v[j]][x] for x in u[i:])
            for w, c in rec(i, j+1).items(): out[(v[j],)+w] += c*wt
        memo[(i, j)] = out; return out
    return rec(0, 0)
def lin(X, Y, C, qq):
    out = defaultdict(lambda: 0)
    for a, ca in X.items():
        for b, cb in Y.items():
            for w, c in shmul(a, b, C, qq).items(): out[w] += ca*cb*c
    return out
def theta_word(w, C, qq):
    X = {(): 1}
    for x in w: X = lin(X, {(x,): 1}, C, qq)
    return X
def adc_power(i, j, n, C, qq):  # ad_c(F_i)^n F_j in free algebra, as dict word->coef
    X = {(j,): 1}
    for _ in range(n):
        deg = defaultdict(int)
        w0 = next(iter(X)); 
        for x in w0: deg[x] += 1
        Y = defaultdict(lambda: 0)
        for w, c in X.items():
            Y[(i,)+w] += c
            Y[w+(i,)] -= c*qq**sum(C[i][x]*m for x, m in deg.items())
        X = Y
    return X
types = {"A1xA1": ([[2, 0], [0, 2]], [1, 1]), "A2": ([[2, -1], [-1, 2]], [1, 1]),
         "B2": ([[2, -2], [-1, 2]], [1, 2]), "G2": ([[2, -3], [-1, 2]], [1, 3])}  # Cartan a_ij, symmetrizer d_i
for name, (A, d) in types.items():
    Cs = [[d[i]*A[i][j] for j in range(2)] for i in range(2)]  # (a_i,a_j) = d_i a_ij
    for i, j in [(0, 1), (1, 0)]:
        n = 1 - A[i][j]
        S = adc_power(i, j, n, Cs, q)
        Th = defaultdict(lambda: 0)
        for w, c in S.items():
            for w2, c2 in theta_word(w, Cs, q).items(): Th[w2] += c*c2
        z = all(sp.simplify(v) == 0 for v in Th.values())
        print(f"{name} Serre ad(F_{i})^{n} F_{j}: Theta(S)=0 symbolically? {z}")
# (2) A3 sorted-word character, exact rationals
from fractions import Fraction as Fr
C = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
def qint_fact(n, z):
    r = Fr(1)
    for k in range(1, n+1): r *= sum(z**e for e in range(k))
    return r
def run(order, qq, Ssym):
    pos = {v: n for n, v in enumerate(order)}
    phi = [[(C[k][l] if pos[k] > pos[l] else 0) + Ssym[k][l] for l in range(3)] for k in range(3)]
    t = [Fr(2), Fr(-5, 3), Fr(7, 4)]
    def Qcorr(nu):  # Q(nu) with Q(a+b)-Q(a)-Q(b) = sigma(a,b), sigma = Ssym
        return sum(Ssym[k][l]*nu[k]*nu[l] for k in range(3) for l in range(k+1, 3)) + sum(Ssym[k][k]*nu[k]*(nu[k]-1)//2 for k in range(3))
    def lam(w):  # sorted-word functional, letters weakly increasing in 'order'
        if any(pos[w[a]] > pos[w[a+1]] for a in range(len(w)-1)): return Fr(0)
        r = Fr(1)
        for k in range(3):
            m = w.count(k); r *= Fr(1)/qint_fact(m, qq**C[k][k])
        return r
    def zeta_free(w):
        r = Fr(1)
        for a in range(len(w)):
            r *= t[w[a]]
            for b in range(a+1, len(w)): r *= qq**phi[w[a]][w[b]]
        return r
    bad = tot = 0
    random.seed(1)
    for _ in range(150):
        w = tuple(random.randrange(3) for _ in range(random.randint(1, 6)))
        nu = [w.count(k) for k in range(3)]
        tw = Fr(1)
        for x in w: tw *= t[x]
        rhs = qq**Qcorr(nu) * tw * sum(Fr(c)*lam(w2) for w2, c in theta_word(w, C, qq).items())
        tot += 1; bad += zeta_free(w) != rhs
    return bad, tot
for qq in [Fr(3, 7), Fr(-5, 2)]:
    for order in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (2, 1, 0)]:
        for Ssym in [[[0]*3]*3, [[1, 2, -1], [2, 0, 3], [-1, 3, -2]]]:
            print("q", qq, "order", order, "sym part", "0" if Ssym[0][0] == 0 else "rand", "zeta_free == q^Q lambda Theta D_t: fails", run(order, qq, Ssym))
# control: lambda_o with the WRONG order phi (reverse) must fail
# NEGATIVE CONTROL: lambda built for order o, but zeta_free uses phi of a different order -> must fail
def run_mismatch(order_lam, order_phi, qq):
    pl = {v: n for n, v in enumerate(order_lam)}; pp = {v: n for n, v in enumerate(order_phi)}
    phi = [[C[k][l] if pp[k] > pp[l] else 0 for l in range(3)] for k in range(3)]
    t = [Fr(2), Fr(-5, 3), Fr(7, 4)]
    bad = 0; random.seed(1)
    for _ in range(150):
        w = tuple(random.randrange(3) for _ in range(random.randint(1, 6)))
        lam = lambda w2: Fr(0) if any(pl[w2[a]] > pl[w2[a+1]] for a in range(len(w2)-1)) else Fr(1)/ (qint_fact(w2.count(0), qq**2)*qint_fact(w2.count(1), qq**2)*qint_fact(w2.count(2), qq**2))
        zf = Fr(1)
        for a in range(len(w)):
            zf *= t[w[a]]
            for b in range(a+1, len(w)): zf *= qq**phi[w[a]][w[b]]
        tw = Fr(1)
        for x in w: tw *= t[x]
        bad += zf != tw*sum(Fr(c)*lam(w2) for w2, c in theta_word(w, C, qq).items())
    return bad
print("CONTROL lam(0,1,2) vs phi(2,1,0): fails", run_mismatch((0, 1, 2), (2, 1, 0), Fr(3, 7)), "/150")
print("CONTROL lam(0,1,2) vs phi(0,2,1): fails", run_mismatch((0, 1, 2), (0, 2, 1), Fr(3, 7)), "/150")
print("CONTROL lam(0,1,2) vs phi(1,0,2) [same orientation of edge 1-2, reversed 0-1]: fails", run_mismatch((0, 1, 2), (1, 0, 2), Fr(3, 7)), "/150")
print("SANITY lam(0,2,1) vs phi(2,0,1) [same orientation: 0 before 1, 2 before 1; 0,2 non-adjacent]: fails", run_mismatch((0, 2, 1), (2, 0, 1), Fr(3, 7)), "/150")

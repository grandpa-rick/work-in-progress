"""Day 206b: machine checks of the four-step proof of W_r.
(A2)  per pair i<j, symmetric F: Y_i Y_j F = t * T_{i-1..1} T_{j-1..2} pi^2 F      [flint exact]
(K)   sigma_m sigma' G = (1+t) sigma^(2) G for G sym in {X1,X2} and in tail           [flint exact]
(K')  sigma^(2) G = sum_{|A|=2} G^(A) prod cross a_ij                                  [exact rational pts]
(E)   final [z^r] assembly (T1+T2)/((1-t)^2(1+t)) == W_r, symbolic e's, r=0..12        [sympy]
(PI)  pi T_k = T_{k+1} pi, k<=m-2                                                      [flint exact]
"""
import sys, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '/home/agent/projects/proofs/scripts/day205')
from k3_fast_pipeline import AHA
import sympy as sp

def Tchain(A, F, hi, lo):  # T_hi ... T_lo F  (T_lo applied first)
    for k in range(lo, hi + 1): F = A.T(F, k)
    return F

ok = True
for m in range(2, 7):
    A = AHA(m); X = A.X
    # random non-symmetric test poly for (PI)
    P = X[0]**2 * X[-1] + 3 * X[0] * X[1 % m] + A.s * X[-1]**3
    for k in range(1, m - 1):
        if A.pi(A.T(P, k)) != A.T(A.pi(P), k + 1): ok = False; print('PI fail', m, k)
    for Fs in [A.e(min(2, m)) * A.e(1), A.e(1)**3 + A.e(m), 1 + 0 * A.s]:
        pi2 = A.pi(A.pi(Fs))
        for i in range(1, m + 1):
            for j in range(i + 1, m + 1):
                lhs = A.Y(A.Y(Fs, j), i)
                rhs = A.t * Tchain(A, Tchain(A, pi2, j - 1, 2), i - 1, 1)
                if lhs != rhs: ok = False; print('A2 fail', m, i, j)
        # (K): sigma_m sigma' G vs (1+t) sigma2 G, G = pi^2 Fs
        G = pi2
        sig2 = 0 * A.s
        for b in range(2, m + 1):
            for a in range(1, b):
                sig2 = sig2 + Tchain(A, Tchain(A, G, b - 1, 2), a - 1, 1)
        sp_ = 0 * A.s
        for b in range(2, m + 1): sp_ = sp_ + Tchain(A, G, b - 1, 2)
        sm = 0 * A.s
        for a in range(1, m + 1): sm = sm + Tchain(A, sp_, a - 1, 1)
        if sm != (1 + A.t) * sig2: ok = False; print('K fail', m)
print('A2 / K / PI :', 'OK' if ok else 'FAIL')

# (E) symbolic assembly
s, t = sp.symbols('s t')
def br(n): return sum(t**k for k in range(n)) if n > 0 else 0
ok2 = True
for r in range(0, 13):
    e = {k: sp.Symbol(f'e{k}') for k in range(0, r + 5)}; e[0] = 1
    E = lambda c, k: c**k * e[k] if k >= 0 else 0
    q1 = (1 - t) * e[1]; q2 = (1 - t) * (e[1]**2 - (1 + t) * e[2])
    # [z^r] of E(cz)/z^j = c^{r+j} e_{r+j}
    T1 = s*(1-t)*( s*(e[1]*q1 - q2)*e[r] + (1-s)*( e[1]*(E(1,r+1)-E(t,r+1)) - (E(t,r+2)-E(1,r+2)+q1*E(1,r+1)) ) )
    T2 = (s-1)*( (s/t)*q1*E(t,r+1) + (1-s/t)*(E(t,r+2)-E(t**2,r+2))/t - s*q1*E(1,r+1) - (1-s)*(E(1,r+2)-E(t,r+2)) )
    tot = sp.simplify((T1 + T2) / ((1-t)**2*(1+t)))
    W = s**2*e[2]*e[r] + s*(1-s)*br(r)*e[1]*e[r+1] + (1-s)*sp.cancel(br(r+2)/(1+t))*(br(r+1) - s*(br(r)-1))*e[r+2]
    if sp.expand(sp.cancel(tot - W)) != 0: ok2 = False; print('E fail r=', r)
print('E (assembly == W_r, r=0..12):', 'OK' if ok2 else 'FAIL')

# (H)-mechanism + (K') at random rational points: sum_{a<b} Phi(Xa,Xb) prod cross == W_r numeric
def evs(vals, k):
    if k < 0 or k > len(vals): return Fr(0)
    return sum((Fr(1) if not c else __import__('math').prod(c)) for c in itertools.combinations(vals, k)) if k else Fr(1)
ok3 = True
random.seed(7)
for m in range(3, 8):
    for r in range(1, 6):
        X = [Fr(v, 7) for v in random.sample(range(2, 200), m)]
        sv = Fr(random.randint(2, 9), 7); tv = Fr(random.randint(2, 9), 5)
        a_ = lambda i, j: (X[i] - tv*X[j]) / (X[i] - X[j])
        K = Fr(0)
        for (a, b) in itertools.combinations(range(m), 2):
            w = Fr(1)
            for i in (a, b):
                for j in range(m):
                    if j not in (a, b): w *= a_(i, j)
            tail = [X[k] for k in range(m) if k not in (a, b)]
            K += w * X[a] * X[b] * evs(tail + [sv*X[a], sv*X[b]], r)
        ee = lambda k: evs(X, k)
        bt = lambda n: sum(tv**k for k in range(n)) if n > 0 else 0
        W = sv**2*ee(2)*ee(r) + sv*(1-sv)*bt(r)*ee(1)*ee(r+1) + (1-sv)*Fr(bt(r+2))/(1+tv)*(bt(r+1) - sv*(bt(r)-1))*ee(r+2)
        if K != W: ok3 = False; print('KW fail', m, r)
print("kernel sum == W_r at random points (m=3..7, r=1..5):", 'OK' if ok3 else 'FAIL')

# (TAU) Lemma 1 assembly: p_2(Y)e_r = e1(Y)(e1(Y)e_r) - 2 e2(Y)e_r, with Thm 3.12, Sub-Lemma Z, W_r
rr = sp.Symbol('r', positive=True, integer=True)
u = t**(rr + 1)
def B(n_shift):  # [r+n_shift] via t^r
    return (1 - t**(rr + n_shift)) / (1 - t)
c211 = s*s**2
c2 = s*(1-s)*s*(1+t) - 2*t*s**2
c11 = (1-s)*B(1)*s + s*(1-s)*(t*B(0)+s) - 2*t*s*(1-s)*B(0)
cr2 = (1-s)**2*B(1)*B(2) + s*(1-s)**2*B(2) - 2*t*(1-s)*B(2)/(1+t)*(B(1) - s*(B(0)-1))
tgt = [s**3, -s**3*(t/s - 1/s + t + 1), s*(1-s**2), -(1-s**2)*B(2)*(u - 1 + s*(1+t))/(1+t)]
ok4 = all(sp.simplify(a - b) == 0 for a, b in zip([c211, c2, c11, cr2], tgt))
print('TAU (Lemma 1 all four coeffs, symbolic r):', 'OK' if ok4 else 'FAIL')

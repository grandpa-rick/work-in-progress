"""Wake 241: general-n bridge from DFK 1704.00154 (negaP)=(6.7) + Lemma dalga (2.16) + Thm genmacthm (2.17):
   M_{N;n}^{-1} M_{N-k;n} = Q^{kn} T^{k(N-k)} M_{k;-n}|_{(Q,T)->(Q^-1,T^-1)}   (printed (1.5), (1.2)); n=-1 is bridge (I).
Checked on monomial symmetric functions of degree <= deg. Control: drop the -n (use M_{k;n})."""
import sympy as sp, itertools
from functools import reduce
Q, T = sp.symbols('Q T')
def prod(l): return reduce(lambda a, b: a*b, l, sp.Integer(1))
def symbasis(x, deg):
    N = len(x); out = []
    for lam in itertools.product(range(deg+1), repeat=N):
        if list(lam) == sorted(lam, reverse=True) and 0 < sum(lam) <= deg:
            out.append(sp.expand(sum(prod([v**e for v, e in zip(p, lam)]) for p in set(itertools.permutations(x)))))
    return out
def M(F, x, a, n, qq, tt):
    N = len(x); r = 0
    for I in itertools.combinations(range(N), a):
        c = prod([(tt*x[i]-x[j])/(x[i]-x[j]) for i in I for j in range(N) if j not in I])
        r += prod([x[i] for i in I])**n*c*F.subs({x[i]: qq*x[i] for i in I}, simultaneous=True)
    return sp.cancel(r)
def MNinv(G, x, n):  # inverse of F -> (x_1..x_N)^n F(Qx)
    N = len(x); return sp.cancel((prod(x)/Q**N)**(-n)*G.subs({v: v/Q for v in x}, simultaneous=True))
for N, deg in [(1, 3), (2, 3), (3, 2)]:
    x = sp.symbols(f'x1:{N+1}'); Fs = [sp.Integer(1)] + symbasis(x, deg)
    for n in (-1, 0, 1):
        for k in range(0, N+1):
            bad = sum(sp.cancel(MNinv(M(F, x, N-k, n, Q, T), x, n) - Q**(k*n)*T**(k*(N-k))*M(F, x, k, -n, 1/Q, 1/T)) != 0 for F in Fs)
            ctl = sum(sp.cancel(MNinv(M(F, x, N-k, n, Q, T), x, n) - Q**(k*n)*T**(k*(N-k))*M(F, x, k, n, 1/Q, 1/T)) != 0 for F in Fs)
            print(f"N={N} n={n:+d} k={k}: identity fails on {bad}/{len(Fs)}; control(M_k;+n) fails on {ctl}/{len(Fs)}")

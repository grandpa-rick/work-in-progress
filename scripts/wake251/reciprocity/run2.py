# Q3/Q4 details: zero structure of F°, real roots of H(j), H(-M/2) at odd M (no lattice count there).
import sys; sys.path.insert(0, '.')
from recip import *
import numpy as np
mmax = int(sys.argv[1]); tot = 0; agg = dict(q1=0, q2e=0, q2o=0, thr=0, M0r=0, roots=0, oc=0, M0le2m=0, intOdd=0, Fodd1=0, Fodd_par1=0)
for m in range(1, mmax + 1):
  for r in range(1, m + 1):
    for kap in itertools.combinations(range(1, m + 1), r):
      if max(kap) != m: continue
      S = setup(kap); P = S[0]; d = len(P); tot += 1; rc = longest_chain(P, boxless)
      Fv = [F(S, N) for N in range(0, 2 * d + 5)]
      Hn = fitpoly([(N, Fv[N]) for N in range(0, 2 * d + 3, 2)])
      assert sp.degree(Hn, x) == d
      agg['q1'] += all(Hn.subs(x, N) == Fv[N] for N in range(len(Fv)))
      sgn = (-1) ** d
      Fo = {M: Fint(S, M) for M in range(1, 2 * d + 5)}
      agg['q2e'] += all(Hn.subs(x, -M) == sgn * Fo[M] for M in Fo if M % 2 == 0)
      # Ehrhart-Macdonald with the ODD constituent (theory says must hold): odd poly through F(odd N)
      Fod = fitpoly([(N, Fv[N]) for N in range(1, 2 * d + 4, 2)])
      agg['q2o'] += all(Fod.subs(x, -M) == sgn * Fo[M] for M in Fo if M % 2 == 1)
      zs = [M for M in Fo if Fo[M] == 0]; M0 = max(zs, default=0)
      agg['thr'] += zs == list(range(1, M0 + 1)); agg['M0r'] += (M0 == rc + 1)
      Hj = sp.Poly(sp.expand(Hn.subs(x, 2 * x)), x)
      # zeros of H at j=-M/2: only even M (j integer) can be, by reciprocity iff M<=M0
      zerosM = [M for M in range(1, 4 * d + 10) if Hj.eval(sp.Rational(-M, 2)) == 0]
      agg['roots'] += zerosM == [M for M in range(2, M0 + 1, 2)]
      oddvals = [Hj.eval(sp.Rational(-M, 2)) for M in range(2 * m + 1, 2 * m + 40, 2)]
      agg['oc'] += all(v != 0 for v in oddvals)
      agg['M0le2m'] += M0 <= 2 * m
      agg['intOdd'] += all(Hj.eval(sp.Rational(-M, 2)).q == 1 for M in range(1, 2 * d + 4, 2))
      rr = sorted(float(z.real) for z in np.roots([float(c) for c in Hj.all_coeffs()]) if abs(z.imag) < 1e-7)
      print(m, kap, 'd', d, 'r', rc, 'M0', M0, 'F°zeros', zs, 'H zeros at -M/2: M=', zerosM,
            'real roots', [round(v, 4) for v in rr], 'H(-M/2) odd M=1,3,5:', [Hj.eval(sp.Rational(-M, 2)) for M in (1, 3, 5)],
            'sign(H(-M/2)) M>=2m+1:', set(np.sign([float(v) * (-1) ** 0 for v in oddvals[:6]])), flush=True)
print('TOTAL', tot, agg)

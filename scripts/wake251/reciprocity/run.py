import sys; sys.path.insert(0, '.')
from recip import *
mmax = int(sys.argv[1])
q1 = q2e = q2o = q3 = q4 = 0; tot = 0; fails = []
rows = []
for m in range(1, mmax + 1):
    for r in range(1, m + 1):
        for kap in itertools.combinations(range(1, m + 1), r):
            if max(kap) != m: continue   # new kappas at this m (kappa subset [m], containing m) -- others seen at smaller m
            S = setup(kap); P = S[0]; d = len(P); tot += 1
            Ns = list(range(0, 2 * d + 5))
            Fv = [F(S, N) for N in Ns]
            Hn = fitpoly([(N, Fv[N]) for N in Ns if N % 2 == 0][: d + 1])  # poly in N; H(j)=Hn(2j)
            ok1 = all(Hn.subs(x, N) == Fv[N] for N in Ns)
            evenok = all(Hn.subs(x, N) == Fv[N] for N in Ns if N % 2 == 0)
            q1 += ok1
            if not ok1: fails.append(('Q1', kap, [(N, Fv[N], Hn.subs(x, N)) for N in Ns if Hn.subs(x, N) != Fv[N]][:3]))
            Ms = list(range(1, 2 * d + 5))
            Fo = [Fint(S, M) for M in Ms]
            sgn = (-1) ** d
            re = all(Hn.subs(x, -M) == sgn * Fo[i] for i, M in enumerate(Ms) if M % 2 == 0)
            ro = all(Hn.subs(x, -M) == sgn * Fo[i] for i, M in enumerate(Ms) if M % 2 == 1)
            q2e += re; q2o += ro
            if not (re and ro): fails.append(('Q2', kap, [(M, Hn.subs(x, -M), sgn * Fo[i]) for i, M in enumerate(Ms) if Hn.subs(x, -M) != sgn * Fo[i]][:4]))
            zeroMs = [M for i, M in enumerate(Ms) if Fo[i] == 0]
            M0 = max(zeroMs, default=0)
            thresh = zeroMs == list(range(1, M0 + 1))
            rc = longest_chain(P, boxless)
            # real roots of H(j) at j = -M/2 for M in 1..2d+4 (scan) vs zeros
            Hj = sp.expand(Hn.subs(x, 2 * x))
            rootsM = [M for M in range(1, 4 * d + 10) if Hj.subs(x, sp.Rational(-M, 2)) == 0]
            ok3 = thresh and rootsM == list(range(1, M0 + 1))
            q3 += ok3
            # Q4: M0 <= 2m and F°(M) > 0 for M >= 2m+1 (within range)
            ok4 = M0 <= 2 * m and all(Fo[i] > 0 for i, M in enumerate(Ms) if M >= 2 * m + 1)
            q4 += ok4
            rows.append((m, kap, d, rc, M0, M0 - rc, thresh, rootsM == list(range(1, M0 + 1)), ok1, re, ro, ok4))
            print(m, kap, 'd', d, 'r', rc, 'M0', M0, 'M0-r', M0 - rc, 'Q1', ok1, 'Q2e/o', re, ro, 'Q3', ok3, 'Q4', ok4, 'Fodd/Feven', Fv[:6], flush=True)
print(f'TOTAL {tot}: Q1 poly {q1}; Q2 even {q2e}, odd {q2o}; Q3 {q3}; Q4 {q4}')
from collections import Counter
print('M0-r distribution', Counter(r[5] for r in rows))
print('fails', fails[:10])

# Wake 252: Crown D test. LG half-integer values H_kappa(a+1/2) vs OG/order-polynomial values Omega(P,a+1).
# H_kappa from the independent Demazure source (m<=4) / file (m=5) via wake251/orderpoly/match.Hs_for.
import sys, itertools
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/reciprocity')
from match import Hs_for, j
from orderpoly import omega_poly, shifted, boxless, t
import recip
import sympy as sp
R = sp.Rational; a = sp.Symbol('a'); AR = range(-6, 7)
mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4

def Om(P):
    return sp.expand(omega_poly(P, boxless).subs(t, a + 1)) if P else sp.Integer(1)

def oddF(kap):  # odd-N constituent of the diag-even quasi-polynomial F, as poly in a with N=2a+1
    S = recip.setup(kap); d = len(S[0])
    return sp.expand(recip.fitpoly([(2*n+1, recip.F(S, 2*n+1)) for n in range(d+1)]).subs(recip.x, 2*a+1))

def const_ratio(A, B):
    if B == 0: return None
    q = sp.cancel(A / B)
    return q if q.is_number else None

stats = {}; def_bump = lambda k: stats.__setitem__(k, stats.get(k, 0) + 1)
for m in range(3, mmax + 1):
    Hs = Hs_for(m); tot = 0
    print(f'===== m={m}')
    for kap, H in sorted(Hs.items(), key=lambda x: (sum(x[0]), x[0])):
        if not kap or max(kap) != m: continue   # genuine LG(m) kappas (smaller max = smaller m)
        tot += 1; lam = sorted(kap, reverse=True); ell = len(lam); d = sum(lam)
        P = shifted(lam); Pnd = [b for b in P if b[0] != b[1]]; Pp = shifted(lam[1:])
        # sanity
        Hh = lambda v: H.subs(j, v)
        Hp = Hs[tuple(sorted(lam[1:]))] if lam[1:] else sp.Integer(1)
        s1 = Hh(R(-1, 2)) == R(1, 2**lam[0]); s2 = Hh(R(-3, 2)) == R(1, 2**lam[0]) * (Hp.subs(j, -1) if lam[1:] else 1)
        def_bump('sanity_-1/2') if s1 else None; def_bump('sanity_-3/2') if s2 else None
        Hhalf = sp.expand(H.subs(j, a + R(1, 2))); Hint = sp.expand(H.subs(j, a))
        Bs = {'Om(P)': Om(P), 'Om(P-diag)': Om(Pnd), "Om(P_kappa')": Om(Pp)}
        Fo = oddF(kap)
        As = {'half': Hhalf, 'INT-control': Hint, 'oddF': Fo}
        for an, A in As.items():
            for bn, B in Bs.items():
                for sc in (1, 2):        # argument scale: B(a+s) and B(2a+s)  [B(x) = Omega(.,x+1)]
                  for s in range(-3, 4):   # shift
                    c = const_ratio(A, B.subs(a, sc * a + s))
                    if c is not None:
                        def_bump(f'{an} = c*{bn}({sc}a{s:+d})'); print(f'  HIT {kap} {an} = {c} * {bn}({sc}a{s:+d})')
        if Fo == Hhalf: def_bump('oddF==Hhalf')
        if m == 3:
            print(f'  kappa={tuple(lam)} l={ell} d={d} sanity({s1},{s2})')
            print('    a      :', list(AR))
            print('    2^l*H(a+1/2):', [2**ell * Hhalf.subs(a, v) for v in AR])
            print('    2^m*H(a+1/2):', [2**m * Hhalf.subs(a, v) for v in AR])
            print('    H(a) [ctrl] :', [Hint.subs(a, v) for v in AR])
            for bn, B in Bs.items(): print(f'    {bn:13s}:', [B.subs(a, v) for v in AR])
            print('    2^l*Om(P)@(2a+1):', [2**0 * Bs['Om(P)'].subs(a, 2*v) for v in AR])
            print('    oddF(2a+1)  :', [Fo.subs(a, v) for v in AR])
    print(f'  m={m}: {tot} kappas')
print('STATS', stats)

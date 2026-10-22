"""Q542 (Clio UID 431) + Q543 support, 2026-10-22. Reads the existing N=4 census coh4_c*.txt (no recompute) and
  (1) counts allowed non-G triples and their P1*-incoherent members at each c;
  (2) checks the non-G data is c-independent (Cor. Z^A_4: singleton support E = {[4] minus ab});
  (3) recomputes the type-C rank-2 census at c = 3..8 directly (triangular solve, no (D)), for Theorem AC transport.
G: consecutive parts differ by >= 2 (theorem-W-sharp L157)."""
import ast, sys
sys.argv = ['x']
from eng import *
from cen_core import analyse

def inG(l): return all(l[i+1] - l[i] >= 2 for i in range(len(l) - 1))
D = {c: {k: d for k, d in (ast.literal_eval(x) for x in open('coh4_c%d.txt' % c) if not x.startswith('#'))}
     for c in range(-2, 9)}
nonG = [k for k in D[6] if not (inG(k[2]) and inG(k[3]))]
for c in D:
    inc = [k for k in D[c] if not D[c][k]['zero'] and not D[c][k]['p1s']]
    print('c=%2d allowed %d nonG %d nonG-incoherent %d G-incoherent %d'
          % (c, len(D[c]), len(nonG), sum(k in nonG for k in inc), sum(k not in nonG for k in inc)))
print('nonG: all nu=[4]:', all(k[4] == (1, 2, 3, 4) for k in nonG),
      '| singleton support = complement of ab:',
      all(EA(4, *k) == [tuple(x for x in range(1, 5) if x not in k[:2])] for k in nonG),
      '| data identical for c=-2..8:',
      all(len({(D[c][k]['P'], D[c][k]['p1s'], D[c][k]['nterms']) for c in D}) == 1 for k in nonG))
print('nonG single-monomial R^A:', sum(D[6][k]['nterms'] == 1 for k in nonG), 'of', len(nonG))
for k in nonG:
    if not D[6][k]['p1s']: print('  incoherent nonG', k, 'P =', D[6][k]['P'])

setup(2)
class SysC(Sys):
    def __init__(s, n):
        s.n = n; g = Fr.C.gens; s.w = [S.const(Fr(g[i])) for i in range(n)]
        s.cache = {}; s.subs = subsets(n); s.down = {}
class R: pass
s = SysC(2)
for c in (3, 4, 5, 6, 8):
    cnt = [0, 0, 0]; bad = []
    for lam in s.subs:
        for mu in s.subs:
            if mu < lam: continue
            for kap in s.subs:
                if not (geq(kap, lam) and geq(kap, mu)): continue
                f = {M: s.wpow(M, c) * s.F(lam, M) * s.F(mu, M) for M in s.subs if geq(kap, M)}
                r = s.coord(f, kap).coef(0); o = R(); o.n = r.n; o.d = r.d; a = analyse(o, 2)
                j = 0 if a['zero'] else (1 if a['p1s'] else 2); cnt[j] += 1
                if j == 2: bad.append((lam, mu, kap, a['P']))
    print('C_2 c=%d zero/coherent/incoherent' % c, cnt)
    if c == 4:
        for b in bad: print('  ', b)

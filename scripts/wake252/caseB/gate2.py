# Wake 252 Case-B gate, corrected.  Two boundary objects on X_k:
#   b_k(j)   = H_k(j)-H_k(j-1) = chi(O_D(j)),  D = div(p_k) = sum m_nu X_nu  (NON-reduced in LG: m_nu=2 off-diagonal)
#   d_k(j)   = chi(O_{dX_k}(j)), dX_k reduced union of codim-1 Schubert = Mayer-Vietoris (IE over covers, meets)
#              = Mobius formula sum_{rho<k} -mu(rho,k) H_rho   [cross-checked IE == MOB]
# True K-class of D: kclass.py (G/B pullback solve).  Term-by-term test uses supp(kclass).
import sys, itertools
from fractions import Fraction as Fr
N = sys.argv[1] if len(sys.argv) > 1 else '5'
sys.argv = sys.argv[:1] + [N]
import gate as G, kclass as KC
def nm(r): return ''.join(map(str, r)) or '0'
Ms = list(range(1, 4 * G.N + 4, 2))
rows = []
for L in range(1, G.N + 1):
  for c in itertools.combinations(range(1, G.N + 1), L):
    k = tuple(sorted(c, reverse=True))
    mob = G.mobius_coeffs(k); ie = G.ie_terms(k); K = KC.kclass(k)
    pts = [Fr(t, 2) for t in range(-4 * G.N - 4, 2 * sum(k) + 4)]
    b = lambda x: G.H(k, x) - G.H(k, x - 1)
    d = lambda x: sum(cc * G.H(r, x) for r, cc in mob.items())
    ieok = all(d(x) == sum(s * G.H(m, x) for s, m in ie) for x in pts)
    kok = all(b(x) == sum(cc * G.H(r, x) for r, cc in K.items()) for x in pts)
    reduced = all(b(x) == d(x) for x in pts)
    a = (G.rp(k) + 1) // 2; top = 2 * a + 1; need = list(range(3, top, 2))
    zH = [M for M in Ms if G.H(k, Fr(-M, 2)) == 0]; zb = [M for M in Ms if b(Fr(-M, 2)) == 0]
    zd = [M for M in Ms if d(Fr(-M, 2)) == 0]
    prop = all(M in zb for M in need); law = zH == list(range(3, top + 1, 2))
    bad = {nm(r): [M for M in need if G.H(r, Fr(-M, 2)) != 0] for r in K}
    bad = {r: v for r, v in bad.items() if v}
    # inductive check: does the odd-zero law for rho (smaller) predict the vanishing?  law-window of rho:
    lawpred = {nm(r): list(range(3, 2 * ((G.rp(r) + 1) // 2) + 2, 2)) for r in K}
    rows.append(dict(k=k, B=G.caseB(k), need=need, prop=prop, law=law, ieok=ieok, kok=kok, red=reduced, tbt=not bad))
    print(k, 'B' if G.caseB(k) else 'A', "r'", G.rp(k), 'zH', zH, 'zb', zb, 'zd(reduced bdry)', zd, 'need', need,
          'prop', prop, 'law', law, 'IE=MOB', ieok, 'Kcls', kok, 'b==d', reduced,
          'K', {nm(r): v for r, v in K.items()}, 'nonvanishing terms', bad, flush=True)
print('IE==MOB', sum(R['ieok'] for R in rows), '/', len(rows), '; Kclass reproduces b', sum(R['kok'] for R in rows),
      '; b == reduced boundary', sum(R['red'] for R in rows))
for lab, f in (('A', lambda R: not R['B']), ('B', lambda R: R['B'])):
    S = [R for R in rows if f(R)]; nt = [R for R in S if R['need']]
    print(f"Case {lab}: {len(S)} kappa, {len(nt)} with nonempty need; propagation {sum(R['prop'] for R in S)}/{len(S)}; "
          f"law {sum(R['law'] for R in S)}/{len(S)}; every K-term vanishes {sum(R['tbt'] for R in nt)}/{len(nt)}; "
          f"b reduced {sum(R['red'] for R in S)}/{len(S)}")

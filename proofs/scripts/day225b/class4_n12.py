# 16 OPEN class-4 pairs n<=12 (classify.py): closed formula vs engine at several numeric t (exact Fractions). grade: computed
import sys; sys.path.insert(0, '../day224')
from fractions import Fraction as Fr
import sympy as sp
from lead2 import second
from hl import parts
from kappa import kappa
from class4 import lead, t
pairs = []
for n in range(4, 13):
    for lam in parts(n):
        if len(lam) != 3: continue
        for mu in parts(n):
            if len(mu) != 2 or mu[0] == mu[1] or kappa(lam, mu) != 1: continue
            x, y = mu
            if y == 1 or x - y == 1 or lam[2] == 1 or lam[0] == x - 1: continue
            pairs.append((lam, mu))
print(len(pairs), 'open pairs', flush=True)
ok = bad = 0
for lam in sorted(set(p[0] for p in pairs), key=sum):
    mus = [p[1] for p in pairs if p[0] == lam]
    closed = {mu: lead(lam, mu) for mu in mus}
    for tv in (Fr(2), Fr(3), Fr(-1, 2)):
        eng = second(*lam, tv)
        for mu in mus:
            g = Fr(str(closed[mu].subs(t, sp.Rational(tv.numerator, tv.denominator)))) == eng.get(mu, 0)
            ok += g; bad += not g
            print(lam, mu, 't=', tv, 'match', g, flush=True)
print('ok', ok, 'bad', bad)

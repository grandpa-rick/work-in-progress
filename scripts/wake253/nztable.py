# Q28 table from out_nz_m*.txt: kappa, r', M*, 2^k1 v (normalised), sign relative to (-1)^|kappa|; integrality of 2^k1 v
import sys, re
from fractions import Fraction as Fr
for m in sys.argv[1:]:
    rows = []
    for L in open(f'out_nz_m{m}.txt'):
        if '|k|' not in L or L.startswith('m='): continue
        f = L.split(); k = f[0]; K = int(f[2]); rp = int(f[8]); Ms = int(f[12]); v = Fr(f[14])
        if rp < 2: continue
        k1 = int(k[0]); nv = v * 2 ** k1
        rows.append((k, rp, Ms, nv, 1 if v * (-1) ** K > 0 else -1))
    print(f'== m={m}: kappa | r\' | M* | 2^k1 * H(-M*/2) | sign vs (-1)^|kappa|')
    for R in rows: print('  ', R[0], R[1], R[2], R[3], '+' if R[4] > 0 else 'EXCEPTION')
    print(f'm={m}: 2^k1 v integral {sum(R[3].denominator == 1 for R in rows)}/{len(rows)}; |2^k1 v|=1 in {sum(abs(R[3]) == 1 for R in rows)}; '
          f'sign exceptions {[R[0] for R in rows if R[4] < 0]}')
# empirical sign rule: sign(v) = (-1)^|kappa| if kappa1 < M*, else (-1)^(|kappa|-kappa1)
print('--- sign rule test: sign = (-1)^|k| if k1<M* else (-1)^(|k|-k1)')
for m in sys.argv[1:]:
    ok = tot = 0; bad = []
    for L in open(f'out_nz_m{m}.txt'):
        if '|k|' not in L or L.startswith('m='): continue
        f = L.split(); k = f[0]; K = int(f[2]); rp = int(f[8]); Ms = int(f[12]); v = Fr(f[14]); k1 = int(k[0])
        if rp < 2: continue
        pred = (-1) ** K if k1 < Ms else (-1) ** (K - k1); tot += 1
        if (v > 0) == (pred > 0): ok += 1
        else: bad.append(k)
    print(f'm={m}: {ok}/{tot}', bad)

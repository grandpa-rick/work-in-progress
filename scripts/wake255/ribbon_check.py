# Wake 255 / Q34: on shapes where T and MN can be compared shape-by-shape (connected theta, odd r), tabulate.
import sys, os
from fractions import Fraction as Fr
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
sys.argv = [sys.argv[0], str(N)]
exec(open(os.path.join(HERE, 'compare_T_MN.py')).read().split('# (1) direct MN')[0])
from qmn import ncomp
tab = Counter()
for k, r in pairs:
    th = boxes(k) - boxes(r); s = len(th)
    if s % 2 == 0: continue
    nd = sum(1 for a, b in th if a == b); cc = ncomp(th)
    t = Tv(k, r); m = MN(r, s).get(k, Fr(0))
    if cc == 1 and nd == 1 and s > 1:
        cls = 'conn, 1 diag box, r>1'
        rel = ('|T|=2^(1-r)|MN|' if abs(t) == abs(m) * Fr(2)**(1-s) else 'magnitude other') + (', same sign' if t*m > 0 else (', opp sign' if t*m < 0 else ', a zero'))
    elif cc == 1 and nd == 0 and s > 1: cls = 'conn, off-diag, r>1'; rel = f'T={t}, MN={m}'
    elif cc == 1: cls = f'conn, {nd} diag'; rel = f'T={t}, MN={m}'
    else: cls = f'{cc} components'; rel = ('T!=0' if t else 'T=0') + ', ' + ('MN!=0' if m else 'MN=0')
    tab[(cls, rel)] += 1
for (c, rl), v in sorted(tab.items()): print(f'{v:4d}  {c:28s} {rl}')

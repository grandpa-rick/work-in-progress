# diagnostics for Q28: per-kappa table (r'>=2): normalised v, sign, candidate counts at L=M*
import sys, itertools
from fractions import Fraction as Fr
n = int(sys.argv[1]); sys.argv = ['nz.py', str(n)]
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): import nz
for R in nz.rows:
    if R['rp'] < 2: continue
    k = R['k']; v = R['v']; L = R['Ms']
    c = R['cand']
    print(nz.nm(k), "r'", R['rp'], 'r', R['r'], 'M*', L, '2^k1*v', v * 2 ** k[0], '2^m*v', v * 2 ** n, 'sgn(-1)^|k|v', 1 if v * (-1) ** sum(k) > 0 else -1,
          '| #strict diag-even P->[M*]', c[('P', 'strict_diag_even')][L], 'diag-odd', c[('P', 'strict_diag_odd')][L],
          'P-row1 strict ->[M*-2]', c[('P-row1', 'strict')][L - 2], 'match(2^k1|v|=diag-even[M*])', abs(v) * 2 ** k[0] == c[('P', 'strict_diag_even')][L])

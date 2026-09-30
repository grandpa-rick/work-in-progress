# exploratory: E_k in rescaled basis b_mu = s^{n(mu)} e_mu. For each (k,mu,nu): val_s [e_nu]E_k e_mu vs n(nu)-n(mu); mod-s matrix entry.
from ds_pointeval import *
from t0_op import op, opcoeffs
from s_valuation import interp
from fractions import Fraction as Fr
import sys
T = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(0)
viol = 0
for n in range(2, int(sys.argv[1])+1):
    for k in range(1, n):
        for mu in parts(n-k):
            D = n*n; sv = [Fr(i+2, 3) for i in range(D+1)]
            data = [opcoeffs(k, mu, x, T) for x in sv]
            row = {}
            for nu in parts(n):
                c = interp(sv, [d.get(nu, Fr(0)) for d in data])
                v = next((i for i, a in enumerate(c) if a != 0), None)
                if v is None: continue
                need = nstat(nu)-nstat(mu)
                if v < need: viol += 1; print('VIOL', k, mu, nu, v, need)
                if v == need: row[nu] = str(c[v])
            print(k, mu, '->', row, flush=True)
print('violations', viol)

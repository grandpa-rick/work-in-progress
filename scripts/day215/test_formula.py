"""Check d_{lam mu}(t) = t^{-n(lam')} sum_nu K_{nu' lam} Ktilde_{nu mu'}(t)  (= t^{-n(lam')} <e_lam, Htilde_{mu'}(x;t)>)
   equivalently d(t) = t^{n(mu')-n(lam')} a_{lam mu'}(1/t), a_{lam rho}(t) = [P_rho(x;t)] e_lam = sum_nu K_{nu' lam} K_{nu rho}(t)."""
import sys
sys.path.insert(0, '/home/agent/projects/scripts/day215')
from analyze import load, conj, parts, nstat
from test_kostka import K, padd, pmul
N = int(sys.argv[1]); ok = tot = 0
for n in range(1, N+1):
    D = load(n)
    for lam in parts(n):
        for mu, d in D[lam].items():
            d = [int(x) for x in d]; s = [0]
            for nu in parts(n): s = padd(s, pmul(K(conj(nu), lam, '1'), K(nu, conj(mu), 'tilde')))
            sh = nstat(conj(lam)); good = all(x == 0 for x in s[:sh]) and s[sh:] == d
            ok += good; tot += 1
            if not good: print('FAIL', lam, mu, d, s)
    print('n', n, 'cumulative', ok, '/', tot, flush=True)

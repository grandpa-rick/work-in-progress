"""P2 test (Wake 246).  For every jump>=2 triple (lam,mu,nu) and every pair a<b (H_ab = {w_a w_b = 1}):
  r(q)   := Res_{H_ab} c^nu_{lam mu}  (full q; local parameter eps, w_b = 1/(w_a(1-eps)), spectators numeric)
  R0     := [q^1] r          (Thm 7 at q=0, proved:  R0 = gamma * Res l0_nu(1_{>=ab} P0),  P0 = phi0_lam phi0_mu)
  Rinf   := [p^1] q^{-e} r   (p = 1/q, e = dim + codim nu - codim lam - codim mu)
  mirror-formula (the q=inf transport of Thm 7):  Rinf =? gamma * Res linf_nu(1_{>=ab} Pinf),  Pinf = phiinf_lam phiinf_mu
mode 'direct' (N<=3): compute r(q) from full Ã and compare Rinf with the formula, R0 with Thm 7.
mode 'formula' (any N): only the two q-free formulas (cheap), to compare the supports Z0 = {R0 = 0}, Zinf = {Rinf = 0}.
Two numeric seeds; a triple counts as 'zero' only if zero at both.  usage: resid_inf.py N direct|formula"""
import sys, random, itertools, time, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake246')
from common import *
N = int(sys.argv[1]); mode = sys.argv[2]; S = subsets(N); w = sp.symbols('w1:%d' % (N+1)); dim = N*(N+1)//2
A = load(N); eps = sp.Symbol('eps'); t0 = time.time()
def point(a, b, seed):
    R = random.Random(seed)
    v = [sp.Rational(R.randint(2, 97), R.randint(2, 97)) + R.randint(0, 3) for _ in range(N)]
    v[b-1] = 1 / (v[a-1] * (1 - eps)); return v
def res(e): return sp.cancel(sp.cancel(e * eps).subs(eps, 0))
triples = [(l, m, nu) for i, l in enumerate(S) for m in S[i:] for nu in S
           if geq(nu, l) and geq(nu, m) and len(nu) - max(len(l), len(m)) >= 2]
out = []; P = lambda s: (print(s, flush=True), out.append(s))
P('N=%d mode=%s: %d jump>=2 triples x %d pairs' % (N, mode, len(triples), N*(N-1)//2))
stats = dict(tot=0, f_inf_ok=0, f_0_ok=0, Z0=0, Zinf=0, both0=0, only0=0, onlyinf=0, prop=0)
detail = []
PAIRS = [tuple(map(int, s.split(","))) for s in sys.argv[3:]] or list(itertools.combinations(range(1, N+1), 2))
for a, b in PAIRS:
    Rv = {}
    for seed in (3, 11):
        v = point(a, b, seed); sub = {w[i]: v[i] for i in range(N)}
        x = v[a-1]; gam = (x - 1/x)**2
        Aq = {k: (A[k].subs(sub, simultaneous=True) if A[k] != 0 else 0) for k in A}
        Ph0 = {k: (sp.cancel(Aq[k].subs(q, 0)) if Aq[k] != 0 else 0) for k in Aq}
        Phi = {k: (lead_inf(Aq[k], dim - sum(k[1])) if Aq[k] != 0 else 0) for k in Aq}
        for l, m, nu in triples:
            e = dim + sum(nu) - sum(l) - sum(m)
            F0 = {M: (Ph0[M, l]*Ph0[M, m] if {a, b} <= set(M) else 0) for M in S}
            Fi = {M: (Phi[M, l]*Phi[M, m] if {a, b} <= set(M) else 0) for M in S}
            form0 = sp.cancel(gam.subs(eps, 0) * res(tri_solve(Ph0, F0, S)[nu]))
            forminf = sp.cancel(gam.subs(eps, 0) * res(tri_solve(Phi, Fi, S)[nu]))
            rec = Rv.setdefault((l, m, nu), {'form0': [], 'forminf': [], 'R0': [], 'Rinf': []})
            rec['form0'].append(form0); rec['forminf'].append(forminf)
        if mode == 'direct':
            for i, l in enumerate(S):
                for m in S[i:]:
                    c = tri_solve(Aq, {M: Aq[M, l]*Aq[M, m] for M in S}, S)
                    for nu in S:
                        if (l, m, nu) not in Rv: continue
                        r = res(c[nu]); e = dim + sum(nu) - sum(l) - sum(m)
                        R0 = sp.cancel(sp.diff(r, q).subs(q, 0))
                        Rinf = sp.cancel(sp.diff(sp.cancel(r.subs(q, 1/p) * p**e), p).subs(p, 0))
                        Rv[l, m, nu]['R0'].append(R0); Rv[l, m, nu]['Rinf'].append(Rinf)
    for tr, rec in Rv.items():
        stats['tot'] += 1
        z0 = all(u == 0 for u in rec['form0']); zi = all(u == 0 for u in rec['forminf'])
        # only count residues on pairs where Thm 6 allows a pole (prop = 'pair is a candidate')
        stats['Z0'] += z0; stats['Zinf'] += zi; stats['both0'] += z0 and zi
        stats['only0'] += z0 and not zi; stats['onlyinf'] += zi and not z0
        if z0 != zi: detail.append(('H%d%d' % (a, b),) + tr + ('R0=0' if z0 else 'Rinf=0',))
        if mode == 'direct':
            stats['f_inf_ok'] += all(sp.cancel(u - v_) == 0 for u, v_ in zip(rec['Rinf'], rec['forminf']))
            stats['f_0_ok'] += all(sp.cancel(u - v_) == 0 for u, v_ in zip(rec['R0'], rec['form0']))
    P('pair (%d,%d) done, %.0fs' % (a, b, time.time() - t0))
P('stats: ' + str(stats))
P('(triple,pair) where exactly one end vanishes: %s' % detail)
open('RESULTS_resid_N%d_%s%s.txt' % (N, mode, ''.join('_H%d%d' % pr for pr in PAIRS) if len(sys.argv) > 3 else ''), 'w').write('\n'.join(out) + '\n')

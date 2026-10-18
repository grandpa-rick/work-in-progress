"""Edge-local version (fast, for N=4): for each moment-graph edge, specialise w to the divisor FIRST:
the free coordinate x = w_a, spectator w's = random rationals (2 seeds), q SYMBOLIC; then build Ã_{M,I}, Ã_{M',I}
(fastA, unpadded) and test that col(M)-col(M') == 0 identically in (x,q).  A fail is certain; a pass is
'generic-spectator' evidence (Schwartz-Zippel style, 2 seeds).  usage: python3 -B psi_gkm_edge.py N"""
import sys, itertools, random, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484'); sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake245')
from engine import subsets, q
from fastA import A as fastA
N = int(sys.argv[1]); S = subsets(N); x = sp.Symbol('x')
psi = lambda v: (1 - q*v**2) / (1 - q)
variants = {'Ahat=Ã*Psi': psi, 'Ã raw': lambda v: 1, 'ctrl Ã*Psi^2': lambda v: psi(v)**2, 'ctrl Ã/Psi': lambda v: 1/psi(v)}
edges = []
for M in S:
    sM = set(M)
    for a in range(1, N+1):
        if a not in sM:
            for sg in (1, -1): edges.append(('2e_a(w=%+d)' % sg, M, tuple(sorted(sM | {a})), {a: sp.Integer(sg)}))
    for a, b in itertools.combinations(range(1, N+1), 2):
        if a not in sM and b not in sM: edges.append(('e_a+e_b', M, tuple(sorted(sM | {a, b})), {a: x, b: 1/x}))
        if a in sM and b not in sM: edges.append(('e_a-e_b', M, tuple(sorted((sM - {a}) | {b})), {a: x, b: x}))
def wvec(fix, seed):
    R = random.Random(seed); w = []
    for k in range(1, N+1):
        w.append(fix[k] if k in fix else sp.Rational(R.randint(2, 60), R.randint(2, 60)) + R.randint(1, 4))
    return w
cnt = {v: {} for v in variants}; fails = {v: [] for v in variants}
for typ, M, M2, fix in edges:
    res = {(v, I): True for v in variants for I in S}
    for seed in (5, 17):
        w = wvec(fix, seed)
        for I in S:
            a1, a2 = fastA(N, M, I, w), fastA(N, M2, I, w)
            for v, f in variants.items():
                d = sp.cancel(sp.together(a1*sp.prod([f(w[m-1]) for m in M]) - a2*sp.prod([f(w[m-1]) for m in M2])))
                if d != 0: res[v, I] = False
    for (v, I), g in res.items():
        c = cnt[v].setdefault(typ, [0, 0]); c[1] += 1; c[0] += g
        if not g: fails[v].append((I, typ, M, M2))
out = []; P = lambda s: (print(s, flush=True), out.append(s))
P('N=%d edge-local: q symbolic, x free, spectators random rational (2 seeds); edges=%d (2e_a split w=+-1), labels=%d' % (N, len(edges), len(S)))
for v in variants:
    P('%-14s ' % v + '  '.join('%s %d/%d' % (k, c[0], c[1]) for k, c in sorted(cnt[v].items())))
    if v.startswith('Ahat'):
        for t in fails[v]: P('    FAIL I=%s %s %s--%s' % t)
open('/home/agent/projects/work-in-progress/scripts/wake245/RESULTS_psi_N%d_edge.txt' % N, 'w').write('\n'.join(out) + '\n')

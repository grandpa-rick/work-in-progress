"""Wake 246: GKM census of transported classes.  For sigma in {bar=(q->1/q, w->1/w), q-only, w-only} and
Psi-choice in {transported Psi^sigma, untransported Psi}, test the LG(N,2N) GKM edge congruences of
   f_I(M) = Ã_{M,I}^sigma * Psi^{tau}(M).
Plus breaking controls: raw Ã^bar (no Psi), Ã^bar * Psi(w;q^2), Ã * Psi(w;q^2).
Edge test: f(M) - f(M') vanishes identically on the divisor (2e_a split into w_a=+1, w_a=-1).
mode 'sym' (N<=3): q and all other w symbolic.  mode 'num' (N=4): q symbolic, free coordinate x symbolic,
spectators random rationals, 2 seeds (a FAIL is certain, a PASS is generic-spectator evidence).
usage: python3 -B bar_gkm_census.py N sym|num"""
import sys, pickle, random, itertools, time, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484')
from engine import q, subsets
N = int(sys.argv[1]); mode = sys.argv[2]; S = subsets(N); w = sp.symbols('w1:%d' % (N+1)); x = sp.Symbol('x')
A = pickle.load(open('/home/agent/projects/work-in-progress/scripts/wake245/A_full_N%d.pkl' % N, 'rb'))
psi = lambda v, qq: (1 - qq*v**2) / (1 - qq)
# transports of Ã, as substitution dicts applied to the stored expression
TR = {'id': {}, 'bar': dict([(q, 1/q)] + [(v, 1/v) for v in w]), 'q-only': {q: 1/q}, 'w-only': {v: 1/v for v in w}}
# variant name -> (transport of Ã, factor chi(v) per m in M)
V = {
 'P1  Ã^bar * Psi^bar      ': ('bar',    lambda v: psi(1/v, 1/q)),
 'P1m Ã^bar * Psi          ': ('bar',    lambda v: psi(v, q)),
 'Cq  Ã(w,1/q)*Psi(w,1/q)  ': ('q-only', lambda v: psi(v, 1/q)),
 'Cqm Ã(w,1/q)*Psi         ': ('q-only', lambda v: psi(v, q)),
 'Cw  Ã(1/w,q)*Psi(1/w,q)  ': ('w-only', lambda v: psi(1/v, q)),
 'Cwm Ã(1/w,q)*Psi         ': ('w-only', lambda v: psi(v, q)),
 'K1  Ã^bar raw (no Psi)   ': ('bar',    lambda v: 1),
 'K2  Ã^bar * Psi(w,q^2)   ': ('bar',    lambda v: psi(v, q**2)),
 'K3  Ã * Psi(w,q^2)       ': ('id',     lambda v: psi(v, q**2)),
 'ref Ã * Psi (Thm Psi)    ': ('id',     lambda v: psi(v, q)),
}
edges = []   # (type, M, M', a, b)
for M in S:
    sM = set(M)
    for a in range(1, N+1):
        if a not in sM:
            for sg in (1, -1): edges.append(('2e_a(w=%+d)' % sg, M, tuple(sorted(sM | {a})), a, sg))
    for a, b in itertools.combinations(range(1, N+1), 2):
        if a not in sM and b not in sM: edges.append(('e_a+e_b', M, tuple(sorted(sM | {a, b})), a, b))
        if a in sM and b not in sM: edges.append(('e_a-e_b', M, tuple(sorted((sM - {a}) | {b})), a, b))
def divisor_point(typ, a, b, seed):
    """values of w on the divisor: free coordinate x; spectators symbolic (sym) or random rationals (num)."""
    R = random.Random(seed); vals = {}
    for k in range(1, N+1):
        vals[k] = w[k-1] if mode == 'sym' else sp.Rational(R.randint(2, 60), R.randint(2, 60)) + R.randint(1, 4)
    if typ.startswith('2e_a'): vals[a] = sp.Integer(b)
    elif typ == 'e_a+e_b': vals[a] = x; vals[b] = 1/x
    else: vals[a] = x; vals[b] = x
    return vals
t0 = time.time()
cache = {}
def Atr(tr, M, I, vals, key):
    k = (tr, M, I, key)
    if k not in cache:
        e = A[M, I]
        if e == 0: cache[k] = sp.Integer(0)
        else:
            e = e.subs(TR[tr], simultaneous=True) if TR[tr] else e
            cache[k] = e.subs({w[i]: vals[i+1] for i in range(N)}, simultaneous=True)
    return cache[k]
cnt = {v: {} for v in V}; fails = {v: [] for v in V}
seeds = (0,) if mode == 'sym' else (5, 17)
for typ, M, M2, a, b in edges:
    res = {(v, I): True for v in V for I in S}
    for sd in seeds:
        vals = divisor_point(typ, a, b, sd); key = (typ, a, b, sd); cache.clear()
        for I in S:
            for v, (tr, chi) in V.items():
                f1 = Atr(tr, M, I, vals, key) * sp.prod([chi(vals[m]) for m in M])
                f2 = Atr(tr, M2, I, vals, key) * sp.prod([chi(vals[m]) for m in M2])
                if sp.cancel(sp.together(f1 - f2)) != 0: res[v, I] = False
    for (v, I), g in res.items():
        c = cnt[v].setdefault(typ, [0, 0]); c[1] += 1; c[0] += g
        if not g: fails[v].append((I, typ, M, M2))
out = []; P = lambda s: (print(s, flush=True), out.append(s))
P('N=%d mode=%s edges=%d labels=%d  (counts = passing (edge,label) pairs / total)  time %.0fs' % (N, mode, len(edges), len(S), time.time()-t0))
for v in V:
    P('%s ' % v + '  '.join('%s %d/%d' % (k, c[0], c[1]) for k, c in sorted(cnt[v].items())))
P('generic-label spot check (fails of K-controls at the largest non-degenerate label):')
for v in V:
    fl = [f for f in fails[v] if 0 < len(f[0]) < N]
    if fl: P('  %s e.g. %s' % (v.strip(), fl[:2]))
open('RESULTS_bar_gkm_N%d_%s.txt' % (N, mode), 'w').write('\n'.join(out) + '\n')

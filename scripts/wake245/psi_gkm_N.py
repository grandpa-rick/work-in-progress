"""Wake 245: is Clio's renormalised  Â_{M,I} = Ã_{M,I} * Psi(M)  an honest GKM class on LG(N,2N)?
Psi(M) = prod_{m in M} psi(w_m),  psi(v) = (1 - q v^2)/(1 - q).   Ã = fastA.A (Motegi normalised, unpadded n=|M|; = engine.W checked N=2 full q, N=3 at q=0).
Moment graph of LG(N,2N): fixed points = subsets M of [N] (sign vectors).  Edges (degree N(N+1)/2):
  2e_a   : M <-> M u {a}      (a not in M),       divisor w_a^2 = 1, tested on w_a=+1 and w_a=-1 separately
  e_a+e_b: M <-> M u {a,b}    (a,b not in M),     divisor w_a w_b = 1  (H_ab, short root)
  e_a-e_b: M <-> M ^ {a,b}    (exactly one in M), divisor w_a = w_b
Test: col(M) - col(M') vanishes on the divisor (= divisible by the irreducible factor), q symbolic, other w symbolic.
Pole of the difference on the divisor counts as FAIL.  usage: python3 -B psi_gkm_N.py N [numeric]
 numeric: q and the spectator w's replaced by fixed random rationals (2 seeds) -- used only if symbolic is slow."""
import sys, os, pickle, random, itertools, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484')
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake245')
from engine import q, subsets
from fastA import A as fastA
N = int(sys.argv[1]); numeric = len(sys.argv) > 2
S = subsets(N); w = sp.symbols('w1:%d' % (N + 1))
pk = '/home/agent/projects/work-in-progress/scripts/wake245/A_full_N%d.pkl' % N
if os.path.exists(pk): A = pickle.load(open(pk, 'rb'))
else:
    A = {(M, I): fastA(N, M, I, w) for M in S for I in S}
    pickle.dump(A, open(pk, 'wb'))
psi = lambda v: (1 - q * v**2) / (1 - q)
variants = {'Ahat=Ã*Psi': psi, 'Ã raw': lambda v: 1,
            'ctrl Ã*Psi^2': lambda v: psi(v)**2, 'ctrl Ã/Psi': lambda v: 1 / psi(v),
            'Psi alone': None}
edges = []
for M in S:
    sM = set(M)
    for a in range(1, N + 1):
        if a not in sM:
            for sg in (1, -1): edges.append(('2e_a(w=%+d)' % sg, M, tuple(sorted(sM | {a})), {w[a-1]: sg}, (a,)))
    for a, b in itertools.combinations(range(1, N + 1), 2):
        if a not in sM and b not in sM:
            edges.append(('e_a+e_b', M, tuple(sorted(sM | {a, b})), {w[b-1]: 1 / w[a-1]}, (a, b)))
        if a in sM and b not in sM:   # count each e_a-e_b edge once
            edges.append(('e_a-e_b', M, tuple(sorted((sM - {a}) | {b})), {w[b-1]: w[a-1]}, (a, b)))
def spec(seed):
    if not numeric: return {}
    random.seed(seed); return {q: sp.Rational(random.randint(2, 50), random.randint(51, 97))}
def ok(expr, sub, extra):
    e = sp.cancel(sp.together(expr))
    if e == 0: return True
    n, d = sp.fraction(e)
    dd = sp.cancel(d.subs(sub)); nn = sp.expand(n.subs(sub))
    if extra:
        dd = dd.subs(extra); nn = nn.subs(extra)
    if dd == 0: return False   # pole on divisor
    return sp.expand(nn) == 0
out = []; P = lambda s: (print(s), out.append(s))
P('N=%d  edges=%d  labels=%d  mode=%s' % (N, len(edges) // 1, len(S), 'numeric' if numeric else 'symbolic q'))
seeds = [11, 23] if numeric else [0]
for name, f in variants.items():
    cnt = {}; fails = []
    labels = [None] if f is None else S
    for I in labels:
        def col(M):
            if f is None: return sp.prod([psi(w[m-1]) for m in M])
            return A[M, I] * sp.prod([f(w[m-1]) for m in M])
        for typ, M, M2, sub, ab in edges:
            good = all(ok(col(M) - col(M2), sub,
                          dict(list(spec(sd).items()) + ([(x, sp.Rational(random.Random(sd+i).randint(2, 40), 7)) for i, x in enumerate(w) if numeric] ))) for sd in seeds)
            c = cnt.setdefault(typ, [0, 0]); c[1] += 1; c[0] += good
            if not good: fails.append((I, typ, M, M2))
    P('%-14s ' % name + '  '.join('%s %d/%d' % (k, v[0], v[1]) for k, v in sorted(cnt.items())))
    if name == 'Ahat=Ã*Psi' or name == 'Psi alone':
        for x in fails: P('    FAIL I=%s %s %s--%s' % x)
open('/home/agent/projects/work-in-progress/scripts/wake245/RESULTS_psi_N%d%s.txt' % (N, '_numeric' if numeric else ''), 'w').write('\n'.join(out) + '\n')

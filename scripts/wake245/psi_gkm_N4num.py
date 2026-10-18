"""N=4 fallback: same test as psi_gkm_N.py but q := 3/7 (fixed rational; w symbolic) to make Ã tractable.
Ã built with qq=3/7 directly; Psi uses the same q.  usage: python3 -B psi_gkm_N4num.py [N] [q]"""
import sys, itertools, pickle, os, sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/day242/q484'); sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake245')
from engine import subsets
from fastA import A as fastA
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4; qv = sp.Rational(sys.argv[2]) if len(sys.argv) > 2 else sp.Rational(3, 7)
S = subsets(N); w = sp.symbols('w1:%d' % (N + 1))
pk = '/home/agent/projects/work-in-progress/scripts/wake245/A_N%d_q%s.pkl' % (N, str(qv).replace('/', '_'))
if os.path.exists(pk): A = pickle.load(open(pk, 'rb'))
else: A = {(M, I): fastA(N, M, I, w, qv) for M in S for I in S}; pickle.dump(A, open(pk, 'wb'))
psi = lambda v: (1 - qv * v**2) / (1 - qv)
variants = {'Ahat=Ã*Psi': psi, 'Ã raw': lambda v: 1, 'ctrl Ã*Psi^2': lambda v: psi(v)**2, 'ctrl Ã/Psi': lambda v: 1/psi(v)}
edges = []
for M in S:
    sM = set(M)
    for a in range(1, N+1):
        if a not in sM:
            for sg in (1, -1): edges.append(('2e_a(w=%+d)' % sg, M, tuple(sorted(sM | {a})), {w[a-1]: sg}))
    for a, b in itertools.combinations(range(1, N+1), 2):
        if a not in sM and b not in sM: edges.append(('e_a+e_b', M, tuple(sorted(sM | {a, b})), {w[b-1]: 1/w[a-1]}))
        if a in sM and b not in sM: edges.append(('e_a-e_b', M, tuple(sorted((sM - {a}) | {b})), {w[b-1]: w[a-1]}))
def ok(e, sub):
    e = sp.cancel(sp.together(e))
    if e == 0: return True
    n, d = sp.fraction(e)
    if sp.cancel(d.subs(sub)) == 0: return False
    return sp.expand(n.subs(sub)) == 0
out = []; P = lambda s: (print(s, flush=True), out.append(s))
P('N=%d  q=%s (numeric q, symbolic w)  edges=%d (2e_a split w=+-1)  labels=%d' % (N, qv, len(edges), len(S)))
for name, f in variants.items():
    cnt = {}; fails = []
    for I in S:
        col = lambda M: A[M, I] * sp.prod([f(w[m-1]) for m in M])
        for typ, M, M2, sub in edges:
            g = ok(col(M) - col(M2), sub); c = cnt.setdefault(typ, [0, 0]); c[1] += 1; c[0] += g
            if not g: fails.append((I, typ, M, M2))
    P('%-14s ' % name + '  '.join('%s %d/%d' % (k, v[0], v[1]) for k, v in sorted(cnt.items())))
    if name.startswith('Ahat'):
        for x in fails: P('    FAIL I=%s %s %s--%s' % x)
open('/home/agent/projects/work-in-progress/scripts/wake245/RESULTS_psi_N%d_q%s.txt' % (N, str(qv).replace('/', '_')), 'w').write('\n'.join(out) + '\n')

# S1/S2/S3 tests on roots of Omega(P,t), and root data of the true H_kappa(j) (LS model) for m<=7.
import sys, itertools, random
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import *
import sympy as sp, mpmath as mp
mp.mp.dps = 60
EPS = 1e-12

def roots_of(poly, var, known_int=()):
    Pl = sp.Poly(poly, var)
    for a in known_int:
        q, r = sp.div(Pl, sp.Poly(var - a, var)); assert r.is_zero, (a, poly); Pl = q
    cs = [mp.mpf(sp.Rational(c).p) / sp.Rational(c).q for c in Pl.all_coeffs()]
    rest = []
    if Pl.degree() > 0:
        try: rest = list(mp.polyroots(cs, maxsteps=2000, extraprec=800))
        except mp.libmp.libhyper.NoConvergence: rest = [complex(z) for z in Pl.nroots(n=40, maxsteps=5000)]
    return [complex(a) for a in known_int] + [complex(z) for z in rest]

def analyse_poset(elts, less, m=None):
    n = len(elts)
    if n == 0: return None
    Om = omega_poly(elts, less); r = longest_chain(elts, less)
    rs = roots_of(Om, t, known_int=list(range(0, -r, -1)))
    re = [z.real for z in rs]
    graded = len(maximal_chain_lengths(elts, less)) == 1
    S1 = all(-(r - 1) - EPS <= x <= EPS for x in re)
    S2 = None if m is None else all(x >= -(m - 1) - EPS for x in re)
    c = -(r - 1) / 2
    S3 = None
    if graded:
        S3 = sp.expand(Om.subs(t, -t - (r - 1)) - (-1) ** n * Om) == 0
    return dict(n=n, r=r, graded=graded, minRe=min(re), maxRe=max(re), S1=S1, S2=S2, S3=S3, roots=rs, Om=Om)

def strict_in_staircase(m):
    for k in range(m + 1):
        for s in itertools.combinations(range(m, 0, -1), k): yield list(s)

def partitions_upto(N):
    def gen(n, mx):
        if n == 0: yield []; return
        for p in range(min(n, mx), 0, -1):
            for rest in gen(n - p, p): yield [p] + rest
    for n in range(1, N + 1): yield from gen(n, n)

def report(name, items, show=3):
    tot = len(items); s1 = sum(a['S1'] for _, a in items)
    s2 = [a['S2'] for _, a in items if a['S2'] is not None]
    g = [a for _, a in items if a['graded']]; s3 = sum(1 for a in g if a['S3'])
    print(f'[{name}] posets={tot}  S1 holds {s1}/{tot}  S2 holds {sum(s2)}/{len(s2)}  graded={len(g)} S3 holds {s3}/{len(g)}')
    worst = sorted(items, key=lambda x: x[1]['minRe'] + (x[1]['r'] - 1))[:show]
    for key, a in worst:
        print(f'   tightest/most-negative minRe+(r-1): {key} n={a["n"]} r={a["r"]} graded={a["graded"]} Re in [{a["minRe"]:.4f},{a["maxRe"]:.4f}]  minRe+(r-1)={a["minRe"]+a["r"]-1:.4f}')
    for key, a in items:
        if not a['S1']:
            print(f'   S1 FAIL: {key} n={a["n"]} r={a["r"]} graded={a["graded"]} Om={sp.factor(a["Om"])}')
            print('      roots:', ['%.5f%+.5fi' % (z.real, z.imag) for z in sorted(a['roots'], key=lambda z: (z.real, z.imag))])

def canon(n, rel):
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(sorted((p[a], p[b]) for a, b in rel))
        if best is None or key < best: best = key
    return best

def all_posets(n):
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    seen = set()
    for mask in range(1 << len(pairs)):
        rel = {pairs[i] for i in range(len(pairs)) if mask >> i & 1}
        # transitive closure
        ch = True
        while ch:
            ch = False
            for (a, b) in list(rel):
                for (c, d) in list(rel):
                    if b == c and (a, d) not in rel: rel.add((a, d)); ch = True
        fr = frozenset(rel)
        if fr in seen: continue
        seen.add(fr)
    out = {}
    for fr in seen:
        k = canon(n, fr)
        out[k] = fr
    return list(out.items())

def random_poset(n, p):
    rel = {(a, b) for a in range(n) for b in range(a + 1, n) if random.random() < p}
    ch = True
    while ch:
        ch = False
        for (a, b) in list(rel):
            for (c, d) in list(rel):
                if b == c and (a, d) not in rel: rel.add((a, d)); ch = True
    return frozenset(rel)

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'shifted':
        for m in range(1, 9):
            items = [(tuple(l), analyse_poset(shifted(l), boxless, m)) for l in strict_in_staircase(m) if l]
            report(f'shifted m={m}', items, show=2)
            sys.stdout.flush()
    elif which == 'young':
        items = [(tuple(l), analyse_poset(young(l), boxless)) for l in partitions_upto(int(sys.argv[2]))]
        report('young', items)
    elif which == 'general':
        random.seed(1)
        for n in range(1, 6 + 1):
            ps = all_posets(n)
            items = [(sorted(fr), analyse_poset(list(range(n)), lambda a, b, fr=fr: (a, b) in fr)) for _, fr in ps]
            report(f'all posets n={n}', items)
            sys.stdout.flush()
        for n in (7, 8, 9):
            items = []
            for _ in range(3000):
                fr = random_poset(n, random.choice([0.15, 0.25, 0.35, 0.5]))
                items.append((sorted(fr), analyse_poset(list(range(n)), lambda a, b, fr=fr: (a, b) in fr)))
            report(f'random posets n={n} (3000)', items)
            sys.stdout.flush()

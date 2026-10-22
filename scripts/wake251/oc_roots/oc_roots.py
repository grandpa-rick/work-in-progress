# Wake 251 compute: root structure of H_kappa(j) = chi(X_kappa, O(j)) on LG(m,2m), kappa subset [m].
# Normalisation: Cor D1 of proofs/2026-10-21-window-is-kodaira-vanishing.md:
#   P^kappa(c) = (-1)^{|kappa|} chi(X_kappa, O(-c/2-1)) = chi(X_kappa, omega_X(c/2+1)),  |kappa| = sum(kappa) = dim X_kappa.
# So a root c0 of P <-> root j0 = -c0/2 - 1 of H; k = c/2+1 = -j is the omega-twist in (OC).
# Odd c <-> half-integer j.  (OC): H has no zero at half-integer j <= -m-1/2.
# Engine: prove249d/dem.py:chi (Demazure character, exact integers) for m<=4;
#         m=5 from prove249d/hp0big2.py (exact localisation, re-run -> hp0big2_n5_rerun.txt), parsed.
import sys, itertools, re
from collections import Counter
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/prove249d')
from dem import chi

j = sp.Symbol('j'); c = sp.Symbol('c')

def H_dem(kap, m):
    d = sum(kap)
    vals = [chi(set(kap), m, x)[0] for x in range(d + 3)]
    H = sp.expand(sp.interpolate(list(zip(range(d + 3), vals)), j))
    assert sp.degree(H, j) == d or d == 0, (kap, H)
    return H

def H_from_P(P, d):
    # P(c) = (-1)^d H(-c/2-1)  =>  H(j) = (-1)^d P(-2j-2)
    return sp.expand((-1) ** d * P.subs(c, -2 * j - 2))

def top_sets(m, shift=0):
    rho = [m - i + shift for i in range(m)]
    longr = Counter(Fr(-r) for r in rho)
    short = Counter(Fr(-(rho[a] + rho[b]), 2) for a in range(m) for b in range(a + 1, m))
    return longr, short

def sub(A, B):
    return all(B[k] >= v for k, v in A.items())

def analyse(H, m, shift=0):
    longr, short = top_sets(m, shift)
    allint = longr + Counter({k: v for k, v in short.items() if k.denominator == 1})
    shorthalf = Counter({k: v for k, v in short.items() if k.denominator == 2})
    if sp.degree(H, j) <= 0:
        return dict(const=True)
    lc, facs = sp.factor_list(H, j)
    rat = Counter(); nonlin = []
    for f, e in facs:
        dg = sp.degree(f, j)
        if dg == 1:
            a, b = sp.Poly(f, j).all_coeffs()
            r = sp.Rational(-b, a); rat[Fr(int(r.p), int(r.q))] += e
        else:
            nonlin.append((f, dg, e))
    ints = Counter({k: v for k, v in rat.items() if k.denominator == 1})
    halfs = Counter({k: v for k, v in rat.items() if k.denominator == 2})
    other = Counter({k: v for k, v in rat.items() if k.denominator > 2})
    nlr = [[complex(z) for z in sp.Poly(f, j).nroots(n=30, maxsteps=2000)] for f, _, _ in nonlin]
    re_ = [float(k) for k in rat] + [z.real for zs in nlr for z in zs]
    d = sp.degree(H, j)
    # symmetry H(j) = +-H(-r-j), r = -2*(sum of roots)/d
    s = sum(sp.Rational(k.numerator, k.denominator) * v for k, v in rat.items())
    r_c = None
    sumroots = -sp.Poly(H, j).all_coeffs()[1] / sp.Poly(H, j).all_coeffs()[0]
    r_c = -2 * sumroots / d
    sym = sp.expand(H.subs(j, -r_c - j) - (-1) ** d * H) == 0
    oc_viol = [k for k in halfs if k <= Fr(-2 * m - 1, 2)]
    return dict(const=False, rat=rat, ints=ints, halfs=halfs, other=other, nonlin=nonlin,
                split=(len(nonlin) == 0),
                H1=sub(halfs, shorthalf), H1strict=sub(halfs, shorthalf) and not other,
                H2=sub(ints, longr), H2all=sub(ints, allint),
                H1set=set(halfs) <= set(shorthalf), H2set=set(ints) <= set(longr),
                minRe=min(re_), maxRe=max(re_),
                strip=all(-(m + 1) < x < 0 for x in re_),
                HPstar=all(x >= -m - 1e-9 for x in re_),
                centre=r_c, sym=sym, oc_viol=oc_viol, nonlin_roots=[sorted(zs, key=lambda z: (z.real, z.imag)) for zs in nlr])

def fmt(cn):
    return '{' + ', '.join(f'{k}' + (f'^{v}' if v > 1 else '') for k, v in sorted(cn.items())) + '}'

def run(m, source):
    print(f'===== m={m}  (variable j: H(j)=chi(X_kappa,O(j)); top-cell index m+1={m+1}) =====')
    subs = [k for r in range(m + 1) for k in itertools.combinations(range(1, m + 1), r)]
    Hs = source(m, subs)
    # sanity: top cell vs Weyl product
    rho = [m - i for i in range(m)]
    W = sp.Integer(1)
    for a in range(m):
        W *= (j + rho[a]) / sp.Integer(rho[a])
        for b in range(a + 1, m):
            W *= (2 * j + rho[a] + rho[b]) / sp.Integer(rho[a] + rho[b])
    top = tuple(range(1, m + 1))
    print('SANITY top cell == Weyl product:', sp.expand(Hs[top] - W) == 0)
    tally = Counter(); res = {}
    for kap in subs:
        H = Hs[kap]; A = analyse(H, m); res[kap] = A
        if A['const']:
            print(kap, 'const', H); tally['const'] += 1; continue
        tally['n'] += 1
        for key in ('split', 'H1', 'H1strict', 'H1set', 'H2', 'H2all', 'H2set', 'strip', 'HPstar', 'sym'):
            tally[key] += A[key]
        tally['halfroot'] += bool(A['halfs']); tally['oc_viol'] += bool(A['oc_viol'])
        nl = '; '.join(f'deg{dg}^{e}: {f}  roots {[f"{z.real:.3f}{z.imag:+.3f}i" for z in zs]}' for (f, dg, e), zs in zip(A['nonlin'], A['nonlin_roots']))
        print(kap, 'dim', sum(kap), '| int', fmt(A['ints']), '| half', fmt(A['halfs']), '| other-rat', fmt(A['other']),
              '| split', A['split'], 'H1', A['H1'], 'H1strict', A['H1strict'], 'H2', A['H2'], 'H2all', A['H2all'], 'H1set', A['H1set'], 'H2set', A['H2set'],
              '| Re in [%.3f, %.3f]' % (A['minRe'], A['maxRe']), 'strip', A['strip'], 'HP*', A['HPstar'],
              '| centre -%s sym %s' % (A['centre'], A['sym']), '| OCviol', A['oc_viol'])
        if nl: print('     nonlinear:', nl)
    print('TALLY m=%d:' % m, dict(tally))
    # ---- controls ----
    # C1 shifted polynomial H(j+1/2): must break H1/H2
    br1 = sum(1 for kap in subs if not res[kap]['const'] and not (lambda A: A['H1'] and A['H2all'])(analyse(sp.expand(Hs[kap].subs(j, j + sp.Rational(1, 2))), m)))
    # C2 reference with shifted rho (rho+1): H2/H1 against wrong top cell
    br2 = sum(1 for kap in subs if not res[kap]['const'] and not (lambda A: A['H1'] and A['H2all'])(analyse(Hs[kap], m, shift=1)))
    # C3 reference of LG(m-1) (wrong ambient)
    br3 = sum(1 for kap in subs if not res[kap]['const'] and not (lambda A: A['H1'] and A['H2all'])(analyse(Hs[kap], m - 1))) if m > 1 else None
    print(f'CONTROLS m={m}: #nonconst={tally["n"]}; (H1&H2all) broken by H(j+1/2): {br1}; by rho+1 reference: {br2}; by LG(m-1) reference: {br3}')
    print(f'  #kappa with half-integer roots (odd-c zeros not from long roots): {tally["halfroot"]}')
    return res, Hs

def src_dem(m, subs):
    return {k: H_dem(k, m) for k in subs}

def src_file(path):
    def f(m, subs):
        out = {}
        for line in open(path):
            mm = re.match(r'^(\(.*?\)) dim (\d+) .* P (.*)$', line.strip())
            if not mm: continue
            kap = eval(mm.group(1)); d = int(mm.group(2))
            P = sp.sympify(mm.group(3), locals={'c': c})
            out[tuple(kap)] = H_from_P(P, d)
        assert set(out) == set(subs), (len(out), len(subs))
        return out
    return f

def crosscheck_hpC3x():
    # independent engine (functional ell_kappa, prove249c/hpCx.py) vs Demazure, lam=mu=empty, n=3
    ok = 0; tot = 0
    for line in open('/home/agent/projects/work-in-progress/scripts/prove249c/hpC3x.txt'):
        mm = re.match(r'^\[\] \[\] (\[.*?\]) d 0 .* P (.*)$', line.strip())
        if not mm: continue
        kap = tuple(eval(mm.group(1))); P = sp.sympify(mm.group(2), locals={'c': c})
        tot += 1; ok += sp.expand(H_from_P(P, sum(kap)) - H_dem(kap, 3)) == 0
    print(f'CROSSCHECK m=3 Demazure vs functional engine (hpC3x.txt, lam=mu=empty): {ok}/{tot}')

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'x': crosscheck_hpC3x()
    elif which in ('3', '4'): run(int(which), src_dem)
    else: run(5, src_file(sys.argv[2]))

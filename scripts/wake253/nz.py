# Wake 253 Q28: (NZ) data for LG(m,2m).  H_kappa from LS-parity model (wake252 gate6.py; verified vs Demazure m<=5, vs BR K-class m<=7).
# For every kappa: odd zero set, W = 2ceil(r'/2)+1 (Thm Z window top), M* = W+2, v = H(-M*/2), middle range odd M in (W, 2m+1) and beyond,
# integer zeros, first integer point past integer window, and candidate counts for normalised |v|.
import sys, itertools, collections
from fractions import Fraction as Fr
n = int(sys.argv[1]); sys.argv = ['gate6.py', str(n)]
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake252/caseB')
import gate6
from gate6 import H, r, rp, nm, shifted, boxless, ideals, longest_chain
def v2(x):   # 2-adic valuation of a nonzero Fraction
    a, b = x.numerator, x.denominator; e = 0
    while a % 2 == 0: a //= 2; e += 1
    while b % 2 == 0: b //= 2; e -= 1
    return e
def layer_counts(P, Lmax, pred):   # # chains of ideals 0=I0<=..<=I_L=P with each layer satisfying pred(layer, v)
    J = ideals(P, boxless); full = frozenset(P)
    sub = {I: [K for K in J if I <= K] for I in J}
    vec = {frozenset(): 1}; out = [1 if not P else 0]
    for v in range(1, Lmax + 1):
        nv = collections.defaultdict(int)
        for I, c in vec.items():
            for K in sub[I]:
                if pred(K - I, v): nv[K] += c
        vec = nv; out.append(vec.get(full, 0))
    return out
def anti(S): return longest_chain(sorted(S), boxless) <= 1 if S else True
def ht2(S): return longest_chain(sorted(S), boxless) <= 2 if S else True
PRED = {
  'strict': lambda S, v: anti(S),
  'strict_diag_even': lambda S, v: anti(S) and (v % 2 == 0 or not any(a == b for a, b in S)),
  'strict_diag_odd': lambda S, v: anti(S) and (v % 2 == 1 or not any(a == b for a, b in S)),
  'fibre_ht2': lambda S, v: ht2(S),
  'weak_diag_odd': lambda S, v: (v % 2 == 1 or not any(a == b for a, b in S)),
}
LMAX = 4 * n + 6
rows = []
for L in range(1, n + 1):
  for c in itertools.combinations(range(1, n + 1), L):
    k = tuple(sorted(c, reverse=True)); rpk = rp(k); W = 2 * ((rpk + 1) // 2) + 1; Ms = W + 2
    oddM = list(range(1, 4 * n + 8, 2)); zH = [M for M in oddM if H(k, Fr(-M, 2)) == 0]
    mid = [M for M in range(W + 2, 2 * n + 1, 2)]; midzero = [M for M in mid if H(k, Fr(-M, 2)) == 0]
    beyond = [M for M in oddM if M > W and H(k, Fr(-M, 2)) == 0]
    v = H(k, Fr(-Ms, 2)); e = -v2(v) if v else None
    zi = [j for j in range(-(sum(k) + 4), 0) if H(k, j) == 0]
    jstar = min(zi) - 1 if zi else -1; hj = H(k, jstar)
    cand = {}
    for nmP, P in (('P', shifted(list(k))), ('P-row1', [(a - 1, b - 1) for (a, b) in shifted(list(k[1:]))])):
        for pn, pr in PRED.items():
            cand[(nmP, pn)] = layer_counts(P, LMAX, pr)
    R = dict(k=k, rp=rpk, r=r(k), W=W, Ms=Ms, v=v, e=e, zH=zH, mid=mid, midzero=midzero, beyond=beyond, zi=zi, jstar=jstar, hj=hj, cand=cand)
    rows.append(R)
    print(nm(k), '|k|', sum(k), 'l', len(k), 'r', R['r'], "r'", rpk, 'W', W, 'M*', Ms, 'v=H(-M*/2)', v,
          '2^e*v', (v * 2 ** e if v else 0), 'e', e, 'sgn', (1 if v > 0 else -1 if v < 0 else 0),
          '2^m v', v * 2 ** n, '| oddzeros', zH, 'midrange', mid, 'midzeros', midzero, 'zeros>W', beyond,
          '| intzeros', zi, 'H(j*)', jstar, hj, flush=True)
rows2 = [R for R in rows if R['rp'] >= 2]
print(f"\nm={n}: kappa with r'>=2: {len(rows2)}; v(M*)!=0 in {sum(R['v'] != 0 for R in rows2)}/{len(rows2)}")
print(f"m={n}: odd zero set exactly [3..W] (all kappa): {sum(R['zH'] == list(range(3, R['W'] + 1, 2)) for R in rows)}/{len(rows)}")
print(f"m={n}: unexpected zeros in middle range (W,2m+1): {[(nm(R['k']), R['midzero']) for R in rows if R['midzero']]}")
print(f"m={n}: any odd zero beyond W up to {4*n+7}: {[(nm(R['k']), R['beyond']) for R in rows if R['beyond']]}")
# sign rule tests
for lab, f in (("(-1)^|k|", lambda R: (-1) ** sum(R['k'])), ("(-1)^(|k|+r'/2...)", lambda R: (-1) ** (sum(R['k']) + (R['W'] - 1) // 2)),
               ("(-1)^((M*-1)/2)", lambda R: (-1) ** ((R['Ms'] - 1) // 2)), ("(-1)^(|k|+l)", lambda R: (-1) ** (sum(R['k']) + len(R['k']))),
               ("(-1)^(|k|+(M*-1)/2)", lambda R: (-1) ** (sum(R['k']) + (R['Ms'] - 1) // 2))):
    print(f"m={n}: sign(v) == {lab}: {sum((R['v'] > 0) == (f(R) > 0) for R in rows2)}/{len(rows2)}")
print(f"m={n}: integer zero set: ", collections.Counter((R['zi'] == list(range(-len(R['zi']), 0)), len(R['zi']) - R['r']) for R in rows).most_common(),
      '(contiguous?, #intzeros - r)')
print(f"m={n}: #int zeros vs candidates: =r' {sum(len(R['zi']) == R['rp'] for R in rows)}, =ceil(r/2) {sum(len(R['zi']) == (R['r']+1)//2 for R in rows)}, "
      f"=floor(r/2) {sum(len(R['zi']) == R['r']//2 for R in rows)}, =ceil(r'/2) {sum(len(R['zi']) == (R['rp']+1)//2 for R in rows)}, =(W-1)/2 {sum(len(R['zi']) == (R['W']-1)//2 for R in rows)}")
# count matching: for each (normalisation, candidate, offset d) with L = M* + d, how many kappa match |norm v| == cand[L]
best = []
for norm in ('2^e', '2^l', '2^m', '2^|k|', '2^k1'):
    for key in rows2[0]['cand']:
        for d in range(-6, 4):
            ok = 0
            for R in rows2:
                v = abs(R['v']); s = {'2^e': 2 ** R['e'] if R['e'] is not None else 1, '2^l': 2 ** len(R['k']), '2^m': 2 ** n,
                                       '2^|k|': 2 ** sum(R['k']), '2^k1': 2 ** R['k'][0]}[norm]
                L = R['Ms'] + d
                ok += 0 <= L < len(R['cand'][key]) and v * s == R['cand'][key][L]
            best.append((ok, norm, key, d))
best.sort(reverse=True)
print(f"m={n}: best uniform count matches (out of {len(rows2)}):", best[:8])
# integer-side reciprocity check: H(-j) == (-1)^|k| #strict g: P->{1..2j-1}, diag even
okr = sum(all(H(R['k'], -j) == (-1) ** sum(R['k']) * R['cand'][('P', 'strict_diag_even')][2 * j - 1] for j in range(1, 2 * n + 2)) for R in rows)
print(f"m={n}: integer reciprocity H(-j) = (-1)^|k| #strict P->[2j-1] diag-even: {okr}/{len(rows)}")

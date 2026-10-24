# Wake 254 Q30 compute.  LG(m,2m): BR K-Pieri coefficients c_rho of [O_{X_kappa}].O^1 (wake252/caseB/brpieri.py, dimension-indexed)
# vs filling counts of theta = kappa/rho (shifted boxes, (row,col), diag row==col).  Also OG checks and H(-j) reciprocity checks.
import sys, itertools, collections
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake252/caseB')
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
import brpieri
from brpieri import boxes, dual, comps, Nprime
from orderpoly import ideals, boxless, shifted
M = int(sys.argv[1])
def strict_parts(n): return [tuple(sorted(c, reverse=True)) for L in range(n + 1) for c in itertools.combinations(range(1, n + 1), L)]
def nm(t): return ''.join(map(str, t)) or '0'
def le(a, b): return a[0] <= b[0] and a[1] <= b[1]
# ---------- fillings ----------
def marked(th, diag_prime_ok):
    """single-letter marked shifted tableaux: labels 0=1', 1=1; rows weakly incr L->R, cols weakly incr T->B,
       1' at most once per row, 1 at most once per column."""
    th = sorted(th); n = 0
    for lab in itertools.product((0, 1), repeat=len(th)):
        L = dict(zip(th, lab)); ok = True
        for (r, c), x in L.items():
            if not diag_prime_ok and r == c and x == 0: ok = False; break
            if (r, c + 1) in L and L[(r, c + 1)] < x: ok = False; break
            if (r + 1, c) in L and L[(r + 1, c)] < x: ok = False; break
        if ok:
            for r in {b[0] for b in th}:
                if sum(1 for (a, c), x in L.items() if a == r and x == 0) > 1: ok = False; break
        if ok:
            for c in {b[1] for b in th}:
                if sum(1 for (a, d), x in L.items() if d == c and x == 1) > 1: ok = False; break
        n += ok
    return n
def strictmaps(th, diag):   # strict order-preserving th -> {0=e<1=o} under P_kappa order restricted; diag in {None,0,1}
    th = sorted(th)
    from orderpoly import longest_chain
    if longest_chain(th, lambda a, b: a != b and le(a, b)) >= 3: return 0
    if len(th) > 14:   # f^-1(e) = down-set I with I, th\I antichains
        anti = lambda S: not any(a != b and le(a, b) for a in S for b in S); n = 0
        for I in ideals(th, lambda a, b: a != b and le(a, b)):
            O = set(th) - I
            if anti(I) and anti(O) and (diag is None or all(b[0] != b[1] for b in (O if diag == 0 else I))): n += 1
        return n
    n = 0
    for lab in itertools.product((0, 1), repeat=len(th)):
        L = dict(zip(th, lab))
        if diag is not None and any(L[b] != diag for b in th if b[0] == b[1]): continue
        if all(L[a] < L[b] for a in th for b in th if a != b and le(a, b)): n += 1
    return n
def ncomp(th): return len(comps(th))
def ndiagcomp(th): return sum(1 for cc in comps(th) if any(r == c for r, c in cc))
def nmax(th): return sum(1 for a in th if not any(a != b and le(a, b) for b in th))
def nmin(th): return sum(1 for a in th if not any(a != b and le(b, a) for b in th))
VARS = {
  'V1 marked P (no diag prime)': lambda th: marked(th, False),
  'V2 marked Q (diag prime ok)': lambda th: marked(th, True),
  'V3e strict e<o, diag=e': lambda th: strictmaps(th, 0),
  'V3o strict e<o, diag=o': lambda th: strictmaps(th, 1),
  'V3* strict e<o, diag free': lambda th: strictmaps(th, None),
  'V4a 2^(comp-diagcomp)': lambda th: 2 ** (ncomp(th) - ndiagcomp(th)),
  'V4b 2^(comp-1)': lambda th: 2 ** (ncomp(th) - 1),
}
def refl(th, n): return {(n - 1 - c, n - 1 - r) for (r, c) in th}   # anti-diagonal reflection of staircase
# ---------- LG ----------
stats = {(v, s): [0, 0, None] for v in VARS for s in ('theta', 'thetaBR')}
supp = {(v, s): [0, 0, None] for v in VARS for s in ('theta', 'thetaBR')}
skipped = collections.Counter(); nterm = 0; reflok = 0; signok = 0; vals = collections.Counter(); KLG = {}
KS = strict_parts(M)
for k in KS:
    K = brpieri.pieri_dim(k, M) if k else {}; KLG[k] = K
    for rho in KS:
        bk, br = set(boxes(k)), set(boxes(rho))
        if not (br < bk): continue
        th = bk - br; thBR = set(boxes(dual(rho, M))) - set(boxes(dual(k, M)))
        c = K.get(rho, 0)
        if c:
            nterm += 1; vals[abs(c)] += 1
            reflok += refl(thBR, M) == th
            signok += (c > 0) == (len(th) % 2 == 1)
        for s, T in (('theta', th), ('thetaBR', thBR)):
            for v, f in VARS.items():
                if not c and len(T) > 14 and not v.startswith('V3'): skipped[v] += 1; continue
                cnt = f(T)
                if c:
                    st = stats[(v, s)]
                    if cnt == abs(c): st[0] += 1
                    else:
                        st[1] += 1
                        if st[2] is None: st[2] = (nm(k), nm(rho), c, cnt, sorted(T))
                # support test: is c==0 iff count==0 (restricted to rim shapes w/ <=2 per row/col for 'count' variants?) -> plain
                sp_ = supp[(v, s)]
                if (cnt != 0) == (c != 0): sp_[0] += 1
                else:
                    sp_[1] += 1
                    if sp_[2] is None: sp_[2] = (nm(k), nm(rho), c, cnt, sorted(T))
print(f'== LG m={M}: {len(KS)} kappa, {nterm} nonzero BR coefficients; |c| histogram {sorted(vals.items())}')
print(f'   theta=kappa/rho equals anti-diagonal reflection of BR nu/lambda: {reflok}/{nterm}; sign = (-1)^(|theta|-1): {signok}/{nterm}')
print('   support-test skips (|theta|>14, off-support, non-V3 variants; per shape-orientation x2):', dict(skipped))
for (v, s), (a, b, ex) in stats.items():
    sa, sb, sex = supp[(v, s)]
    print(f'   {v:30s} on {s:8s}: |c|==count {a}/{a+b}  first fail {ex} | support (c!=0 <=> count!=0 over all rho<kappa) {sa}/{sa+sb} first fail {sex}')
# ---------- OG ----------
sys.argv = ['og.py', str(M), 'none']
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake253')
import og
KOG = {}; t1 = t2 = t3 = t4 = 0; ex = []
for k in KS:
    if not k: KOG[k] = {}; continue
    K = og.kog_pieri_dim(k, M, 1); KOG[k] = K
    P = set(boxes(k)); cor = {a for a in P if not any(a != b and le(a, b) for b in P)}
    thetas = []
    for rho, c in K.items():
        th = P - set(boxes(rho)); thetas.append(frozenset(th))
        t1 += abs(c) == 1; t2 += (c > 0) == (len(th) % 2 == 1); t4 += 1
    want = {frozenset(S) for L in range(1, len(cor) + 1) for S in itertools.combinations(sorted(cor), L)}
    ok = set(thetas) == want and len(thetas) == len(want); t3 += ok
    if not ok and len(ex) < 3: ex.append((nm(k), sorted(cor), [sorted(t) for t in thetas]))
nk = len(KS) - 1
print(f'== OG m={M}: {t4} nonzero KOG coefficients: |c|=1 {t1}/{t4}; sign (-1)^(|theta|-1) {t2}/{t4}; '
      f'supp thetas == nonempty subsets of corners of P_kappa: {t3}/{nk}  fails {ex}')
# ---------- H(-j) via recursion (*) ----------
def Hrec(Kd, J):   # integer values H_kappa(j), j in [-J, J], from H(0)=1 and H(j)-H(j-1)=sum c_rho H_rho(j); kappa in size order
    Hv = {}
    for k in sorted(KS, key=sum):
        h = {0: 1}
        b = lambda j: sum(c * Hv[r][j] for r, c in Kd[k].items())
        for j in range(1, J + 1): h[j] = h[j - 1] + b(j)
        for j in range(0, -J - 1, -1): h[j - 1] = h[j] - b(j)
        Hv[k] = h
    return Hv
def layer_count(k, nvals, pred):   # # maps P_kappa -> {1..nvals}, layers f^-1(v) antichains (strict) with pred(layer,v)
    P = shifted(list(k)); J = ideals(P, boxless); full = frozenset(P)
    anti = lambda S: not any(a != b and le(a, b) for a in S for b in S)
    vec = {frozenset(): 1}
    for v in range(1, nvals + 1):
        nv = collections.defaultdict(int)
        for I, c in vec.items():
            for Kk in J:
                if I <= Kk and anti(Kk - I) and pred(Kk - I, v): nv[Kk] += c
        vec = nv
    return vec.get(full, 0) if P else 1
JJ = M + 2
HL = Hrec(KLG, JJ); HO = Hrec(KOG, JJ)
# cross-check positive side against existing models (LS-parity for LG, Omega for OG)
sys.argv = ['gate6.py', str(M)]; import gate6
posL = all(HL[k][j] == gate6.H(k, j) for k in KS for j in range(0, 4)) ; posO = all(HO[k][j] == og.omega_vals(k, max(sum(k), 3))[j] for k in KS if k for j in range(0, 4))
print(f'== recursion (*) H values agree with LS-parity model (LG) j=0..3: {posL}; with Omega model (OG) j=0..3: {posO}')
okL = 0
for k in KS:
    good = all(HL[k][-j] == (-1) ** sum(k) * layer_count(k, 2 * j - 1, lambda S, v: v % 2 == 0 or not any(a[0] == a[1] for a in S))
               for j in range(1, JJ + 1))
    okL += good
    if not good: print('LG H(-j) FAIL', nm(k), [HL[k][-j] for j in range(1, JJ + 1)])
print(f'== LG: H(-j) == (-1)^|k| #strict maps P_k->[2j-1] even on diagonal, j=1..{JJ}: {okL}/{len(KS)}')
for d in (-1, 0, 1, 2):
    okO = sum(all(HO[k][-j] == (-1) ** sum(k) * (layer_count(k, j - 1 + d, lambda S, v: True) if j - 1 + d >= 0 else (1 if not k else 0))
                  for j in range(1, JJ + 1)) for k in KS)
    print(f'== OG: H(-j) == (-1)^|k| Omega_bar(P_k, j-1{d:+d}) [strict maps to [j-1{d:+d}]], j=1..{JJ}: {okO}/{len(KS)}')

# Wake 257 Q41: is the OG hyperplane identity  H(j)-H(j-1) = sum_rho c_rho H_rho(j)  blind to swapping
# kappa <-> kappa' when P_kappa ~= P_kappa' as UNMARKED posets (e.g. (3) <-> (2,1))?
# Coefficients c_rho = coeff of [O_{X_rho}] in [O_{X_kappa}].O^1 on OG(m,2m+1), computed FOUR ways:
#   BRh : BR eq.(hsumB) literal  B(theta,1)=sum_phi (-1)^{|th|-|phi|} h(N^-(phi), d(phi)-1), theta = rho^vee/kappa^vee (Cor 4.4)
#   KOG : BR KOG-tableau rule (wake253/og.py kog_pieri_dim)
#   TRUE: K-class via pull-back to G/B + type-B Demazure characters (wake253/og.py true_kclass; independent of BR) [m<=TRUEMAX]
#   LEM : draft Lemma og-corners: c_rho=(-1)^{|th|-1} iff th=D(kappa)\D(rho) nonempty set of maximal elements
# H^OG_kappa(j) = Omega(P_kappa, j+1) (og.py omega model; cross-checked vs type-B Demazure chi_dem at small j here).
# Usage: python3 q41.py M TRUEMAX
import sys, itertools, collections, random, time
M = int(sys.argv[1]); TRUEMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 4
sys.argv = [sys.argv[0], str(M)]
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake253')
import og
from fractions import Fraction as Fr
from math import comb
from og import strict_parts, boxes, dual, nm, H, chi_dem, omega_vals, kog_pieri_dim, true_kclass, is_rim
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake252/caseB')
import brpieri  # LG coefficients (for the negative control)

def h(a, b): return sum((-1) ** j * 2 ** (a - j) * comb(a, j) for j in range(0, min(a, b) + 1)) if b >= 0 else 0
def se_corners(th): return {(r, c) for (r, c) in th if (r + 1, c) not in th and (r, c + 1) not in th}
def se_rim(th): return {(r, c) for (r, c) in th if not any(r2 > r and c2 > c for (r2, c2) in th)}
def ncomp(th):
    th = set(th); n = 0
    while th:
        st = [th.pop()]; n += 1
        while st:
            r, c = st.pop()
            for nb in ((r + 1, c), (r - 1, c), (r, c + 1), (r, c - 1)):
                if nb in th: th.remove(nb); st.append(nb)
    return n
def Bh(th, p):
    cor = se_corners(th); base = th - cor; tot = 0
    for L in range(len(cor) + 1):
        for S in itertools.combinations(sorted(cor), L):
            phi = base | set(S)
            tot += (-1) ** (len(th) - len(phi)) * h(max(ncomp(phi) - 1, 0), len(se_rim(phi)) - p)
    return tot
def brh_pieri_dim(k, n):
    lam = dual(k, n); bl = boxes(lam); out = {}
    for nu in strict_parts(n):
        bn = boxes(nu)
        if bl < bn:
            v = Bh(bn - bl, 1)
            if v: out[dual(nu, n)] = v
    return out
def le(a, b): return a[0] <= b[0] and a[1] <= b[1]
def lem_pieri(k, n):
    D = boxes(k); mx = [x for x in D if not any(y != x and le(x, y) for y in D)]; out = {}
    for L in range(1, len(mx) + 1):
        for A in itertools.combinations(mx, L):
            rest = D - set(A); rho = tuple(sorted(collections.Counter(r for r, c in rest).values(), reverse=True))
            assert boxes(rho) == rest
            out[rho] = (-1) ** (L - 1)
    return out
# ---------- unmarked poset isomorphisms ----------
def isos(k1, k2, limit=10 ** 6):
    A, B = sorted(boxes(k1)), sorted(boxes(k2))
    if len(A) != len(B): return []
    def inv(S, x): return (sum(le(y, x) for y in S), sum(le(x, y) for y in S))
    iA = {x: inv(A, x) for x in A}; iB = {y: inv(B, y) for y in B}
    if sorted(iA.values()) != sorted(iB.values()): return []
    A.sort(key=lambda x: (x[0] + x[1], x)); out = []
    def bt(i, mp, used):
        if len(out) >= limit: return
        if i == len(A): out.append(dict(mp)); return
        x = A[i]
        for y in B:
            if y in used or iB[y] != iA[x]: continue
            if all(le(x, z) == le(y, mp[z]) and le(z, x) == le(mp[z], y) for z in mp):
                mp[x] = y; used.add(y); bt(i + 1, mp, used); del mp[x]; used.discard(y)
    bt(0, {}, set()); return out
def marked_ok(phi): return all((x[0] == x[1]) == (y[0] == y[1]) for x, y in phi.items())
def shape(S): return tuple(sorted(collections.Counter(r for r, c in S).values(), reverse=True))
def transport(exp, k, phi):  # rho -> rho' with D(rho') = phi(D(rho))
    out = {}
    for rho, c in exp.items():
        img = {phi[x] for x in boxes(rho)}; r2 = shape(img); assert boxes(r2) == img
        out[r2] = c
    return out

t0 = time.time(); KS = [k for k in strict_parts(M) if k]
print(f'===== m={M} =====')
# coefficients, 3 or 4 ways
EXP = {}; agree = collections.Counter()
LAMS = []
random.seed(1)
for _ in range(6 * 2 ** M): par = random.randint(0, 1); LAMS.append(tuple(2 * random.randint(-og.RR, og.RR) + par for _ in range(M)))
cache = {}
for k in KS:
    e1, e2, e4 = brh_pieri_dim(k, M), kog_pieri_dim(k, M, 1), lem_pieri(k, M)
    agree['BRh==KOG'] += e1 == e2; agree['BRh==LEM'] += e1 == e4
    if M <= TRUEMAX:
        e3 = true_kclass(k, M, LAMS, cache); agree['BRh==TRUE'] += e1 == e3
    EXP[k] = e1
print('coefficient agreement over', len(KS), 'kappa:', dict(agree))
# H model vs Demazure for the shapes of interest
JS = [Fr(t, 2) for t in range(-16, 25)]
def ident(k, exp, Hk):  # H_k(j)-H_k(j-1) == sum c H_rho(j) on grid JS
    return all(Hk(x) - Hk(x - 1) == sum(c * H(q, x) for q, c in exp.items()) for x in JS)
okI = sum(ident(k, EXP[k], lambda x, k=k: H(k, x)) for k in KS)
print(f'identity (BRh coeffs x Omega model) holds on half-int grid j in [-8,12]: {okI}/{len(KS)}')
# (a)+(b) for (3),(2,1)
if M >= 3:
    a, b = (3,), (2, 1)
    for k in (a, b):
        print(f'(a) kappa={nm(k)}: expansion', {nm(q): v for q, v in EXP[k].items()},
              ' Demazure chi j=0..3:', [chi_dem(k, M, j) for j in range(4)], ' Omega:', list(omega_vals(k, 3)[:4]))
    Ph = isos(a, b); print('   unmarked isos (3)->(2,1):', len(Ph), ' marked-compatible:', sum(marked_ok(p) for p in Ph))
    print('   H_(3)==H_(2,1) as polynomials (grid):', all(H(a, x) == H(b, x) for x in JS))
    print('(b) literal swap  H_(3) with (2,1)-expansion:', ident(a, EXP[b], lambda x: H(a, x)),
          '| H_(2,1) with (3)-expansion:', ident(b, EXP[a], lambda x: H(b, x)))
    for p in Ph:
        T = transport(EXP[a], a, p)
        print('    transported (3)-expansion via iso:', {nm(q): v for q, v in T.items()}, ' == (2,1)-expansion:', T == EXP[b],
              '| relabelled swap identity holds:', ident(b, T, lambda x: H(b, x)))
# (c) all unmarked-isomorphic pairs
pairs = []; bysize = collections.defaultdict(list)
for k in KS: bysize[sum(k)].append(k)
for s, L in bysize.items():
    for k1, k2 in itertools.combinations(L, 2):
        Ph = isos(k1, k2)
        if Ph: pairs.append((k1, k2, Ph))
print(f'(c) unmarked-isomorphic pairs kappa!=kappa\' at m={M}: {len(pairs)}')
for k1, k2, Ph in pairs:
    Heq = all(H(k1, x) == H(k2, x) for x in JS)
    lit = EXP[k1] == EXP[k2]
    tr = [transport(EXP[k1], k1, p) == EXP[k2] for p in Ph]
    sw = ident(k1, EXP[k2], lambda x: H(k1, x)) and ident(k2, EXP[k1], lambda x: H(k2, x))
    print(f'   {nm(k1)} ~ {nm(k2)}: #isos {len(Ph)} marked-compatible {sum(marked_ok(p) for p in Ph)}; H equal {Heq}; '
          f'expansions literally equal {lit}; transported==target for all isos {all(tr)} (any {any(tr)}); literal swap identity {sw}; '
          f'exp1 {dict((nm(q), v) for q, v in EXP[k1].items())} exp2 {dict((nm(q), v) for q, v in EXP[k2].items())}')
# ---------- negative controls ----------
print('--- negative controls ---')
if M >= 3:
    a, b = (3,), (2, 1)
    bad = dict(EXP[b]); q0 = next(iter(bad)); bad[q0] = -bad[q0]
    print('NC1 sign-flipped (2,1)-expansion fed to identity for (2,1): holds =', ident(b, bad, lambda x: H(b, x)), '(expect False)')
    La, Lb = brpieri.pieri_dim(a, M), brpieri.pieri_dim(b, M)
    print('NC2 LG coefficients (BR Cor 5.4) (3):', {nm(q): v for q, v in La.items()}, ' (2,1):', {nm(q): v for q, v in Lb.items()})
    PhL = isos(a, b)
    print('    LG transported (3)-exp == (2,1)-exp for some iso:', any(transport(La, a, p) == Lb for p in PhL), '(expect False: comparator can detect a difference)')
    print('    LG coeffs fed into OG identity for (3): holds =', ident(a, La, lambda x: H(a, x)), '(expect False)')
# random non-iso bijection control for the transport/iso checker
k1 = max(KS, key=sum); bx = sorted(boxes(k1)); perm = bx[1:] + bx[:1]
phi = dict(zip(bx, perm)); isiso = all(le(x, y) == le(phi[x], phi[y]) for x in bx for y in bx)
print('NC3 cyclic-shift bijection of', nm(k1), 'accepted as iso by order test:', isiso, '(expect False); isos() self-count:', len(isos(k1, k1)))
print(f'time {time.time()-t0:.1f}s')

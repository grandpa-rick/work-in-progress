# Check, over ALL shifted skew nu/lam in rho_n (not only those hit by q30), that
#  (A) BR Cor 5.4 h-sum C(theta,1)  [proved in BR]
#  (B) (-1)^{|theta|-1} #KLG-tableaux content {1}, BR Def before Cor 5.6, conditions (i)-(v) implemented LITERALLY
#  (C) Rick Q30 closed form on reflected theta: 0 unless no 3-chain & no diag box strictly BELOW another box (reflected:
#      diag box strictly ABOVE... handled by reflecting), else 2^{#isolated off-diagonal boxes}
# agree.  Also checks BR Thm 5.5 recursion at p=1 directly.
import sys, itertools
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake252/caseB')
from brpieri import boxes, C, se_rim
n = int(sys.argv[1])
def le(a, b): return a[0] <= b[0] and a[1] <= b[1]
def lt(a, b): return a != b and le(a, b)
def klg(th):
    if se_rim(th) != th: return 0
    th = sorted(th); cnt = 0
    for lab in itertools.product((0, 1), repeat=len(th)):  # 0 = 1', 1 = 1
        L = dict(zip(th, lab)); ok = True
        for a in th:
            for b in th:
                if a == b: continue
                if a[0] == b[0] and a[1] < b[1] and not L[a] < L[b]: ok = False   # (i) rows strict
                if a[1] == b[1] and a[0] < b[0] and not L[a] < L[b]: ok = False   # (ii) cols strict
                sw = b[0] >= a[0] and b[1] <= a[1]                               # b weakly SW of a
                if sw and L[a] == 1 and not L[a] >= L[b]: ok = False             # (iii)
                if sw and L[a] == 0 and not L[a] <= L[b]: ok = False             # (iv)
            if a[0] == a[1] and L[a] == 0: ok = False                            # (v)
        cnt += ok
    return cnt
def rick(th):   # th in kappa/rho orientation; diag forced e = bottom
    th = list(th)
    if any(lt(a, b) and lt(b, c) for a in th for b in th for c in th): return 0
    if any(d[0] == d[1] and lt(x, d) for d in th for x in th): return 0
    iso = [a for a in th if a[0] != a[1] and not any(lt(a, b) or lt(b, a) for b in th)]
    return 2 ** len(iso)
def refl(th): return {(n - 1 - c, n - 1 - r) for (r, c) in th}
parts = [tuple(sorted(c, reverse=True)) for L in range(n + 1) for c in itertools.combinations(range(1, n + 1), L)]
tot = okAB = okAC = 0; bad = []
for nu in parts:
    for lam in parts:
        bn, bl = boxes(nu), boxes(lam)
        if not (bl < bn): continue
        th = bn - bl; tot += 1
        A = C(th, 1); B = (-1) ** (len(th) - 1) * klg(th); Cc = (-1) ** (len(th) - 1) * rick(refl(th))
        okAB += A == B; okAC += A == Cc
        if (A != B or A != Cc) and len(bad) < 5: bad.append((nu, lam, A, B, Cc))
print(f'n={n}: {tot} skew shapes; h-sum == KLG count {okAB}/{tot}; h-sum == Rick closed form {okAC}/{tot}; bad {bad}')

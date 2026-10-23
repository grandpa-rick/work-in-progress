# Support pattern of [O_D] = [O_{X_kappa}].O^1 (dimension-indexed), from computed classes (gate outputs).
# theta = kappa/rho as shifted skew diagram.  Candidate rule (R): theta nonempty, no two boxes in strict SE relation (rim),
# every row and every column of theta has <= 2 boxes; coefficient = (-1)^{|theta|-1} * #labelings with {1',1}:
# rows/cols strictly increasing, diagonal boxes constrained (variant 'unp': unprimed, 'pr': primed).
import re, ast, sys, itertools, collections
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
def r(k): return longest_chain(shifted(list(k)), boxless) if k else 0
def rp(k): return r(k[1:])
def bx(l): return {(i, i + j) for i, li in enumerate(l) for j in range(li)}
def leq(a, b): return len(a) <= len(b) and all(a[i] <= b[i] for i in range(len(a)))
def labelings(th, diag):
    th = sorted(th); n = 0
    for lab in itertools.product((0, 1), repeat=len(th)):     # 0 = 1', 1 = 1
        L = dict(zip(th, lab))
        if any((r_, c + 1) in L and not (L[(r_, c)] == 0 and L[(r_, c + 1)] == 1) for (r_, c) in L): continue
        if any((r_ + 1, c) in L and not (L[(r_, c)] == 0 and L[(r_ + 1, c)] == 1) for (r_, c) in L): continue
        if any(r_ == c and L[(r_, c)] != (1 if diag == 'unp' else 0) for (r_, c) in L): continue
        n += 1
    return n
def rule(th, diag):
    if not th: return 0
    if any(r2 > r1 and c2 > c1 for (r1, c1) in th for (r2, c2) in th): return 0
    return (-1) ** (len(th) - 1) * labelings(th, diag)
def parse(f):
    for L in open(f):
        if L.startswith('('):
            k = ast.literal_eval(re.match(r'(\(.*?\))', L).group(1))
            K = ast.literal_eval(re.search(r"'?K'? (\{.*?\})", L).group(1))
            yield k, {(tuple(int(ch) for ch in s) if s != '0' else ()): v for s, v in K.items()}
if __name__ == '__main__':
    n = int(sys.argv[2]); data = list(parse(sys.argv[1]))
    allstrict = [tuple(sorted(c, reverse=True)) for L in range(n + 1) for c in itertools.combinations(range(1, n + 1), L)]
    for diag in ('unp', 'pr'):
        ok = sum(all(K.get(rho, 0) == rule(bx(k) - bx(rho), diag) for rho in allstrict if leq(rho, k) and rho != k) for k, K in data)
        print(f'rule (R) diag={diag}: exact K-class {ok}/{len(data)}')
    # statistics
    depth = collections.Counter(); rowmax = colmax = 0; coeffs = collections.Counter(); drops = []
    for k, K in data:
        for rho, v in K.items():
            th = bx(k) - bx(rho); depth[len(th)] += 1; coeffs[v] += 1
            rows = collections.Counter(a for a, _ in th); cols = collections.Counter(b for _, b in th)
            rowmax = max(rowmax, max(rows.values())); colmax = max(colmax, max(cols.values()))
            if len(k) >= 2 and rp(rho) - rp(k) <= -2: drops.append((k, rho, rp(rho) - rp(k), v))
    print('|kappa/rho| distribution:', sorted(depth.items()), '| max boxes per row', rowmax, 'per col', colmax, '| coeffs', sorted(coeffs.items()))
    print("min r'(rho)-r'(kappa) over supp, l(kappa)>=2:", sorted(collections.Counter(
        min(rp(rho) - rp(k) for rho in K) for k, K in data if len(k) >= 2).items()))
    print("pairs with r'(rho)-r'(kappa) <= -2 (kappa, rho, drop, coeff):")
    for d in drops: print('  ', d, 'theta rows', sorted(collections.Counter(a for a, _ in bx(d[0]) - bx(d[1])).items()))

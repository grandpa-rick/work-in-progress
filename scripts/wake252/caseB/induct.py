# From out_gate2_m5.txt: for each kappa, does the odd-zero LAW of every rho in supp[O_D] (inductive hypothesis, rho smaller)
# cover the needed points M in [3, 2ceil(r'(k)/2)-1]?  Window(rho) = [3, 2ceil(r'(rho)/2)+1].  Also: does supp contain rho with l(rho)<=1?
import re, ast, sys
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
def rp(k): return longest_chain(shifted(list(k[1:])), boxless) if len(k) > 1 else 0
tot = dict(A=[0, 0], B=[0, 0]); fails = []
for L in open(sys.argv[1] if len(sys.argv) > 1 else 'out_gate2_m5.txt'):
    if not L.startswith('('): continue
    k = ast.literal_eval(re.match(r'(\(.*?\))', L).group(1)); case = L.split()[len(str(k).split())]
    K = ast.literal_eval(re.search(r"'K' ?(\{.*?\})", L).group(1)) if "'K'" in L else ast.literal_eval(re.search(r"K (\{.*?\})", L).group(1))
    need = list(range(3, 2 * ((rp(k) + 1) // 2) + 1, 2))
    if not need: continue
    ok = True
    for s in K:
        r = tuple(int(ch) for ch in s) if s != '0' else ()
        top = 2 * ((rp(r) + 1) // 2) + 1
        if not all(3 <= M <= top for M in need): ok = False; fails.append((k, s, need, top))
    tot[case][0] += ok; tot[case][1] += 1
print('law-of-rho covers need (term-by-term induction closes):', {c: f'{a}/{b}' for c, (a, b) in tot.items()})
for f in fails: print('FAIL', f)
# support statistics for the inductive lemma
import collections
st = collections.Counter(); mins = []
for L in open(sys.argv[1] if len(sys.argv) > 1 else 'out_gate2_m5.txt'):
    if not L.startswith('('): continue
    k = ast.literal_eval(re.match(r'(\(.*?\))', L).group(1)); K = ast.literal_eval(re.search(r"K (\{.*?\})", L).group(1))
    rs = [tuple(int(ch) for ch in s) if s != '0' else () for s in K]
    if len(k) >= 2:
        mins.append((min(rp(r) for r in rs) - rp(k), max(sum(k) - sum(r) for r in rs), min(len(r) for r in rs) >= 2, k))
print("min over supp of r'(rho)-r'(k):", sorted(collections.Counter(m[0] for m in mins).items()),
      '| max box-depth of supp:', max(m[1] for m in mins), '| all supp have l>=2 (l(k)>=2):', all(m[2] for m in mins))
print('kappa with r\'(rho)-r\'(k) = -2 somewhere:', [m[3] for m in mins if m[0] <= -2])

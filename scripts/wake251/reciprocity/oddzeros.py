# Test: odd-M zeros of H(-M/2) are exactly M=3,5,..,2*ceil(r'/2)+1, r' = longest chain of shifted(kappa minus its largest part).
import re, ast, sys
sys.path.insert(0, '/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
ok = tot = 0
for L in open('out_run2_m6.txt'):
    if L.startswith('TOTAL'): continue
    kap = ast.literal_eval(re.search(r'(\(.*?\))', L).group(1)); zs = ast.literal_eval(re.search(r'M= (\[.*?\])', L).group(1))
    odd = [z for z in zs if z % 2]
    rest = sorted(kap, reverse=True)[1:]; rp = longest_chain(shifted(rest), boxless) if rest else 0
    pred = list(range(3, 2 * ((rp + 1) // 2) + 2, 2)); tot += 1; ok += odd == pred
    if odd != pred: print('MISS', kap, odd, pred)
print(f'odd-zero law: {ok}/{tot}')

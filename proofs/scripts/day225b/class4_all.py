# compare closed class-4 leads vs engine data (newcases_n10.log) for ALL kappa=1 two-part pairs n<=10, all 3 orderings. grade: computed
import sys, re; sys.path.insert(0, '../day224')
import sympy as sp
from class4 import lead, t
ok = bad = 0
for line in open('../day224/newcases_n10.log'):
    m = re.match(r'(\([\d, ]*\)) -> (\([\d, ]*\)) : (.*?)  \|', line)
    if not m: continue
    lam, mu, p = eval(m.group(1)), eval(m.group(2)), sp.sympify(m.group(3)).subs(sp.Symbol('t'), t)
    res = [sp.expand(lead(lam, mu, i) - p) == 0 for i in range(3)]
    g = all(res); ok += g; bad += not g
    print(lam, mu, res)
print('ok', ok, 'bad', bad)

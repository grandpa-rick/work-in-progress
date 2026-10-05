from closed import *
import re
ok = bad = 0
for line in open('newcases_n10.log'):
    m = re.match(r'(\([\d, ]*\)) -> (\([\d, ]*\)) : (.*?)  \|', line)
    if not m: continue
    lam, mu, p = eval(m.group(1)), eval(m.group(2)), sympy.sympify(m.group(3))
    if mu[1] != 1: continue
    c = cN1(lam); e = sympy.symbols('e')
    ser = sympy.expand(sympy.cancel(c.subs(s, 1+e)))
    c0, c1, c2 = [sympy.factor(sympy.cancel(ser.coeff(e, i))) for i in range(3)]
    g = c0 == 0 and c1 == 0 and sympy.expand(c2 - p) == 0
    ok += g; bad += not g
    print(lam, mu, 'v>=2:', c0 == 0 and c1 == 0, ' lead match:', sympy.expand(c2-p) == 0)
print('ok', ok, 'bad', bad)

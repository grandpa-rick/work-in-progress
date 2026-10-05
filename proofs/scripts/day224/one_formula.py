"""Lead_{(b,c,1),(x,y)} closed formula (a=1 ordering): Gamma_1 diagonal-only + Thm 7.1. kappa=1, x>y>=1, x+y=b+c+1."""
import sympy, re
t = sympy.symbols('t')
def q(m): return sum(t**i for i in range(m))
def L(a, b): return -q(a+b)*q(a*b)/(q(a)*q(b))
def lead1(b, c, x, y):
    n1 = b+c; m, M = min(b, c), max(b, c); S = {x, y}
    tot = 0
    for i in range(1, m):
        tot -= L(b-i, c-i) * (q(n1+1-i)*int(i in S) + q(i+1)*int((i+1) in S))
    for j in range(1, m+1): tot -= q(j+1)*int((j+1) in S)
    for j in range(M+1, n1): tot += q(j+1)*int((j+1) in S)
    return sympy.factor(sympy.cancel(tot))
ok = bad = 0
for f in ['newcases_n10.log', 'table3_n8.log']:
    for line in open(f):
        m = re.match(r'(\([\d, ]*\)) -> (\([\d, ]*\)) : (.*?)  \|', line)
        if not m: continue
        lam, mu, p = eval(m.group(1)), eval(m.group(2)), sympy.sympify(m.group(3))
        if lam[2] != 1 or len(mu) != 2: continue
        v = lead1(lam[0], lam[1], *mu); g = sympy.expand(v - p) == 0; ok += g; bad += not g
        print(lam, mu, g, v)
print(ok, bad)
def lead1_simpl(b, c, x, y):
    m, M = min(b, c), max(b, c); I = lambda c_: 1 if c_ else 0
    v = q(x)*(I(x >= M+2) - I(x <= m+1) - (L(b-y, c-y) if y < m else 0)) \
        - q(y)*(I(2 <= y <= m+1) + (L(b-y+1, c-y+1) if 2 <= y <= m else 0))
    return sympy.factor(sympy.cancel(v))
ok2 = bad2 = 0
from hl import parts
from kappa import kappa
for n in range(4, 16):
    for lam in parts(n):
        if len(lam) != 3 or lam[2] != 1: continue
        for mu in parts(n):
            if len(mu) != 2 or mu[0] == mu[1] or kappa(lam, mu) != 1: continue
            g = sympy.expand(lead1(lam[0], lam[1], *mu) - lead1_simpl(lam[0], lam[1], *mu)) == 0
            ok2 += g; bad2 += not g
            if not g: print('SIMPL FAIL', lam, mu)
print('simplified == sum form:', ok2, bad2)

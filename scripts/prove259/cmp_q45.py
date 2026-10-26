# compare combinatorial P (verify.py) with Wake 259 linear-solve P = S C^{-1} S (out_q45_m3.txt), all entries
import re, sympy as sp, sys
sys.argv=['verify.py','3','2']
exec(open('verify.py').read().split('def mul')[0])
X=sp.symbols('x1:4')
comb={(k,r):sp.expand(sum(v*sp.Mul(*[X[t]**e[t] for t in range(3)]) for e,v in c.items())) for (k,r),c in P.items() if c}
old={}
for L in open('../wake259/out_q45_m3.txt'):
    mm=re.match(r'\s+P\[(.*?)\] = (.*?)   \(at',L)
    if mm:
        a,b=mm.group(1).split('),(') if '),(' in mm.group(1) else None
        key=(eval(a+')'),eval('('+b)); old[key]=sp.sympify(mm.group(2))
print(len(comb),len(old), all(sp.simplify(comb.get(k,0)-old.get(k,0))==0 for k in set(comb)|set(old)))

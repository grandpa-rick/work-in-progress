import sys, sympy as sp, re
c=sp.Symbol('c')
for fn,n in [('../prove249c/hpC3x.txt',3),('../prove249c/hpC4x_0.txt',4),('../prove249c/hpC4x_1.txt',4)]:
    tot=ok=okint=dz=0; bad=[]
    for line in open(fn):
        m=re.match(r'(\[.*?\]) (\[.*?\]) (\[.*?\]) d (\d+) deg (\d+).* P (.*)$',line.strip())
        if not m: continue
        lam,mu,kap=[eval(m.group(i)) for i in (1,2,3)]; d=int(m.group(4)); P=sp.sympify(m.group(6))
        tot+=1
        if d>0: continue
        dz+=1
        U=2*n-2; v=P.subs(c,U+2); vo=P.subs(c,U+1)
        s=(-1)**(sum(lam)+sum(mu))
        if sp.sign(v)==s and sp.sign(vo)==s: ok+=1
        else: bad.append((lam,mu,kap))
        if all(P.subs(c,2*k).is_integer for k in range(-3,12)): okint+=1
    print(fn,'total',tot,'d=0',dz,'sign (-1)^{|lam|+|mu|} at U+1,U+2:',ok,'integer at even c:',okint,'bad',bad[:5])

import sys, re, sympy as sp
sys.path.insert(0,'../prove247')
import mp
c_=sp.Symbol('c'); N=int(sys.argv[1]); cs=eval(sys.argv[2])
inst=[]
for fn in sys.argv[3:]:
  for line in open(fn):
    m=re.match(r'(\[.*?\]) (\[.*?\]) (\[.*?\]) d',line)
    if not m: continue
    inst.append((tuple(tuple(eval(m.group(i))) for i in (1,2,3)),sp.sympify(line.split(' P ',1)[1])))
s1=mp.Sys(N,21); s2=mp.Sys(N,22)
for cc in cs:
    zX=set(); zP=set()
    for (lam,mu,kap),P in inst:
        if P.subs(c_,cc)==0: zP.add((lam,mu,kap))
        if mp.C(s1,lam,mu,kap,cc)==0 and mp.C(s2,lam,mu,kap,cc)==0: zX.add((lam,mu,kap))
    print('c',cc,'#X=0',len(zX),'#P=0',len(zP),'X\\P',len(zX-zP),'P\\X',len(zP-zX),flush=True)

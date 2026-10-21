import sys, re, sympy as sp
c=sp.Symbol('c'); t=sp.Symbol('t')
U=int(sys.argv[1]); bad=0; tot=0; chkbad=0; zeros=0
for fn in sys.argv[2:]:
  for line in open(fn):
    if 'ZERO' in line: zeros+=1; print(line.strip()); continue
    if ' P ' not in line: continue
    tot+=1
    if 'chk False' in line: chkbad+=1; print('CHK',line.strip())
    P=sp.sympify(line.split(' P ',1)[1])
    Q=sp.Poly(sp.expand(P.subs(c,U+1+t)),t).all_coeffs()
    nz=[x for x in Q if x!=0]
    if not (Q[-1]!=0 and (all(x>0 for x in nz) or all(x<0 for x in nz))): bad+=1; print('DESC FAIL',line.strip())
print('total',tot,'desc fail',bad,'chk fail',chkbad,'zero-lines',zeros)

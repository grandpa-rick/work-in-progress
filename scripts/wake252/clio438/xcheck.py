# cross-check own one-head Psi1 against Rick's engine q532.Psi (agreement of definitions)
import sys,random,sympy as sp; sys.path.insert(0,'../../prove249')
import contextlib,io
with contextlib.redirect_stdout(io.StringIO()): import q532
from tie23_defs import Psi1
random.seed(7); n=0
for S in [(2,3,4),(1,2,3,4),(1,5,7)]:
  for s in (1,2,3):
    for t in (-2,0,1):
      x=q532.pt(8); x={k:sp.Rational(v.numerator,v.denominator) for k,v in x.items()}
      n+= Psi1(s,t,S,x)==q532.Psi((s,),(t,),S,x)
print('agree',n,'/27')

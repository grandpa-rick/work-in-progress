import sys; sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/day242/q484'); sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/day244')
from engine import *
from compare import koncki, lab, S, t, s, y
from gkm import gkm, pts
w1,w2=sp.symbols('w1 w2'); d={w1:1/s,w2:1/t}
psi=lambda v:(1-q*v**2)/(1-q)
def Psi(M):
  r=1
  for m in M: r*=psi([w1,w2][m-1])
  return r
print('S',S)
for I in S:
  col={lab[M]:(W('C',2,I,M)[0]*Psi(M)).subs(d).subs(q,-y) for M in S}
  print('Psi*Motegi',lab[I],gkm(col))
print('Psi alone',gkm({lab[M]:sp.sympify(Psi(M)).subs(d).subs(q,-y) for M in S}))
lam=sp.Rational(1,4); c=sp.ceiling(lam); c2=sp.ceiling(2*lam); al=2*c-c2
A=sp.Matrix([[ (W('C',2,I,M)[0]*Psi(M)).subs(d).subs(q,-y) for I in S] for M in S])
A0=sp.Matrix([[ W('C',2,I,M)[0].subs(d).subs(q,-y) for I in S] for M in S])
K=sp.Matrix([[ sp.sympify(koncki(lab[I],lab[M],c,c2,al)) for I in S] for M in S])
for name,B in [('Motegi',A0),('Psi*Motegi',A)]:
  C=(B.inv()*K).applyfunc(lambda e: sp.factor(sp.simplify(e)))
  print('C:',name,'-> Koncki (rows J, cols I, order',[lab[I] for I in S],')')
  for i in range(4): print('   ',[C[i,j] for j in range(4)])
C=(A0.inv()*A).applyfunc(lambda e: sp.factor(sp.simplify(e.subs(y,-q))))
print('C: Motegi -> Psi*Motegi (y=-q, s,t)')
for i in range(4): print('   ',[C[i,j] for j in range(4)])

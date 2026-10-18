import sys; sys.path.insert(0,'.'); sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/day242/q484')
from cf import *; from engine import subsets, geq
N=4; w=list(sp.symbols('w1:5')); x,eps=sp.symbols('x eps')
a,b=1,2; nu=(1,2,4)
down=[M for M in subsets(N) if geq(nu,M)]
order=sorted(down,key=lambda T:(len(T),sum(T)))
for wt in (False,True):
  sub={w[a-1]:x, w[b-1]:1/(x*(1-eps))}
  A={(M,I):sp.cancel(F(I,M,w,wt).subs(sub)) for M in down for I in down}
  # l_nu(e_M) via triangular solve on down-set; residue
  def coords(f):
    c={}
    for M in order:
      v=f.get(M,0)-sum(c[I]*A[M,I] for I in c)
      c[M]=sp.cancel(v/A[M,M])
    return c
  for M in down:
    if a in M and b in M:
      r=sp.factor(sp.cancel((coords({M:1})[nu]*eps)).subs(eps,0))
      print('wt',wt,'rho',M,r)
  for Mpp in [(3,),(4,)]:
      print('wt',wt,'P',Mpp, sp.factor(F((3,),Mpp,w,wt).subs(sub).subs(eps,0)))

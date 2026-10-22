import sys, itertools, sympy as sp
sys.path.insert(0,'/home/agent/projects/work-in-progress/scripts/wake251/orderpoly')
from orderpoly import shifted, boxless, longest_chain
K=lambda tau,t:(tau-t)/(tau*t-1)
def calH(lam,c,T):
    l=len(lam); t=sp.symbols('t1:%d'%(l+1)); expr=1
    for s in range(1,l+1):
        ts=t[s-1]
        rho=ts**(2*c-1-s)*sp.Mul(*[K(tau,ts) for tau in T])*sp.Mul(*[K(t[r-1],ts) for r in range(1,s)])
        expr*=ts*(1-rho)/((ts**2-1)*(ts-1)**lam[s-1])
    for s in range(l,0,-1): expr=sp.residue(sp.together(expr),t[s-1],1)
    return sp.simplify(expr)
taus=sp.symbols('tau1:4')
if __name__=="__main__":
   for kap in [(1,),(2,),(3,),(4,),(1,2),(1,3),(2,3),(1,4),(2,4),(3,4),(1,2,3)]:
      lam=sorted(kap,reverse=True); r=longest_chain(shifted(lam),boxless); rest=lam[1:]; rp=longest_chain(shifted(rest),boxless) if rest else 0
      out=[]
      for q in range(0,3):
          zs=[c for c in [sp.Rational(n,2) for n in range(1,13)] if calH(lam,c,taus[:q])==0]
          out.append((q,[str(z) for z in zs]))
      print(kap,'r',r,"r'",rp, out, flush=True)

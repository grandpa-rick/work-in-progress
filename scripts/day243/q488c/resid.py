import sympy as sp,pickle,sys,itertools
N=int(sys.argv[1]); w=sp.symbols('w1:%d'%(N+1)); x=sp.Symbol('x')
res=pickle.load(open('res_N%d.pkl'%N,'rb'))
for (l,m,nu),(j,c0,c1,lp0,lp1) in res.items():
  if j<2: continue
  a,b=sorted(nu)[:2]
  h=1-1/(w[a-1]*w[b-1])
  r=sp.factor(sp.cancel(sp.cancel(c1*h).subs(w[b-1],1/w[a-1])))
  g=(w[a-1]-1/w[a-1])**2
  print(l,m,nu,'H',(a,b),' Res/gamma =',sp.factor(sp.cancel(r/g)))

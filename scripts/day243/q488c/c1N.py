import sympy as sp,sys,time,pickle,itertools,os
from engine import subsets, geq
from engine3 import A01
N=int(sys.argv[1]); S=subsets(N); w=sp.symbols('w1:%d'%(N+1)); t=time.time()
fn='A01_N%d.pkl'%N
if os.path.exists(fn): A0,A1=pickle.load(open(fn,'rb'))
else:
  A0={};A1={}
  for M in S:
    for I in S: A0[M,I],A1[M,I]=A01(N,M,I)
  pickle.dump((A0,A1),open(fn,'wb'))
print('A built',time.time()-t,flush=True)
order=sorted(S,key=lambda T:(len(T),sum(T)))
def solve(vec):
  c={}
  for nu in order:
    rhs=vec[nu]-sum(c[k]*A0[nu,k] for k in c if A0[nu,k]!=0 and c[k]!=0)
    c[nu]=sp.factor(sp.cancel(rhs/A0[nu,nu]))
  return c
def longpoles(e):
  if e==0: return set()
  den=sp.fraction(sp.cancel(e))[1]; out=set(); other=[]
  for f,m in sp.factor_list(den)[1]:
    if not f.free_symbols: continue
    if sp.Poly(f,*w).is_monomial: continue
    hit=False
    for a,b in itertools.combinations(range(1,N+1),2):
      if sp.cancel(f/(w[a-1]*w[b-1]-1)).is_number: out.add((a,b,m)); hit=True
    if not hit: other.append(sp.factor(f))
  return out, other
res={}
for l in S:
  for m in S:
    if order.index(m)<order.index(l): continue
    c0=solve({M:A0[M,l]*A0[M,m] for M in S})
    Phi={M:A0[M,l]*A1[M,m]+A1[M,l]*A0[M,m]-sum(c0[k]*A1[M,k] for k in S if c0[k]!=0) for M in S}
    c1=solve(Phi)
    for nu in S:
      if not (geq(nu,l) and geq(nu,m)):
        assert c0[nu]==0 and c1[nu]==0, ('support',l,m,nu); continue
      j=len(nu)-max(len(l),len(m))
      lp0=longpoles(c0[nu]); lp1=longpoles(c1[nu])
      res[l,m,nu]=(j,c0[nu],c1[nu],lp0,lp1)
      print(l,m,nu,'jump',j,'c0=0' if c0[nu]==0 else 'c0!=0','poles c0',lp0,'poles c1',lp1,flush=True)
pickle.dump(res,open('res_N%d.pkl'%N,'wb'))
print('time',time.time()-t)

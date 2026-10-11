# Prop 2.1 (Subset formula) independent control: realise Cherednik-type operators Y_i directly from
# Demazure-Lusztig T_i and the PRINTED pi F = x_1 F(x_2,...,x_N, s x_1), and compare t^{-C(k,2)} e_k(Y) F with the
# printed subset formula E_k F on symmetric F, N=3. Hikita's exact conventions are not re-derived here, so we search
# a small family of standard conventions and report those for which the Y_i commute AND reproduce E_k for k=1,2,3.
import itertools, sympy as sp
s,t=sp.symbols('s t'); N=3; x=sp.symbols('x1:%d'%(N+1))
def swap(f,i): return f.subs({x[i]:x[i+1],x[i+1]:x[i]},simultaneous=True)
def T_a(f,i): return sp.cancel(t*f+(t*x[i]-x[i+1])/(x[i]-x[i+1])*(swap(f,i)-f))
def T_b(f,i): return sp.cancel(t*f+(t*x[i+1]-x[i])/(x[i+1]-x[i])*(swap(f,i)-f))
def pi_printed(f,sv=s):
    return sp.expand(x[0]*f.subs({x[j]:(x[j+1] if j<N-1 else sv*x[0]) for j in range(N)},simultaneous=True))
def esym(vals,k): return sum(sp.Mul(*c) for c in itertools.combinations(vals,k))
def E_subset(k,f):
    tot=0
    for A in itertools.combinations(range(N),k):
        cA=sp.Mul(*[(x[i]-t*x[j])/(x[i]-x[j]) for i in A for j in range(N) if j not in A])
        XA=sp.Mul(*[x[i] for i in A])
        tot+=cA*XA*f.subs({x[i]:s*x[i] for i in A},simultaneous=True)
    return sp.cancel(tot)
tests=[esym(x,1),esym(x,2),esym(x,1)**2]
Etab={(k,i):E_subset(k,f) for k in (1,2,3) for i,f in enumerate(tests)}
def run(Tf,pos,seq,e,rec,pi):
    Tinv=lambda f,i: sp.cancel((Tf(f,i)-(t-1)*f)/t)
    Tp=Tf if e==1 else Tinv
    order=list(range(N-1)) if seq=='asc' else list(range(N-1))[::-1]
    def Y1(f):
        g=f
        if pos=='piLeft':      # Y1 = pi o T_{order[0]} o T_{order[1]} ...  (rightmost acts first)
            for j in reversed(order): g=Tp(g,j)
            return pi(g)
        else:
            g=pi(g)
            for j in reversed(order): g=Tp(g,j)
            return g
    Ys=[Y1]
    for i in range(N-1):
        prev=Ys[-1]
        if rec==0: Ys.append(lambda f,i=i,prev=prev: sp.cancel(t*Tinv(prev(Tinv(f,i)),i)))
        else:      Ys.append(lambda f,i=i,prev=prev: sp.cancel(Tf(prev(Tf(f,i)),i)/t))
    g=x[0]**2*x[1]
    comm=all(sp.cancel(Ys[i](Ys[j](g))-Ys[j](Ys[i](g)))==0 for i in range(N) for j in range(i+1,N))
    if not comm: return comm,None
    res=[]
    for k in (1,2,3):
        ok=True
        for ii,f in enumerate(tests):
            tot=0
            for S in itertools.combinations(range(N),k):
                g2=f
                for i in S: g2=Ys[i](g2)
                tot+=g2
            if sp.cancel(tot/t**sp.binomial(k,2)-Etab[(k,ii)])!=0: ok=False
        res.append(ok)
    return comm,res
hits=0
for Tn,Tf in [('T_a',T_a),('T_b',T_b)]:
  for pos in ('piLeft','piRight'):
    for seq in ('asc','desc'):
      for e in (1,-1):
        for rec in (0,1):
          comm,res=run(Tf,pos,seq,e,rec,pi_printed)
          if comm: print(Tn,pos,seq,'e=',e,'rec=',rec,'commute=True','match k=1,2,3:',res)
          if res and all(res):
              hits+=1
              # negative control: same convention with pi using s->1 (no q-shift) must FAIL to match
              _,res2=run(Tf,pos,seq,e,rec,lambda f: pi_printed(f,1))
              print('   negative control (pi without s):', res2)
print('conventions reproducing the printed subset formula:',hits)

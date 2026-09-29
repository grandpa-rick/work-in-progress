"""Day 209: check the closed-form step map (a)+(b) of the proof against
   day208 gf_recursion.step output (pickled Gamma_1..3), and against (TC) weights.
Step map on a term W e_b E(g)E(d), g=t^i z, d=t^j w, output times (1-t):
 (a) e_{b+1}E(g)E(d):  W s^2 t^{-i-j} (1-t^{b+1})
 (b) n<=b: e_n * W * [ -kg g^{n-b-1}(E(tg)E(d) - t^n E(g)E(d)) - kd d^{n-b-1}(E(g)E(td) - t^n E(g)E(d)) ]
 kg=(g-sz)(g-sw)/(g(g-d)), kd=(d-sz)(d-sw)/(d(d-g)).  Then divide by (1-t)[k]."""
import sympy as sp, pickle, sys
s,t,z,w=sp.symbols('s t z w')
def step(G,k):
    new={}
    def add(key,v): new[key]=new.get(key,0)+v
    for (b,i,j),Wt in G.items():
        g,d=t**i*z,t**j*w
        kg=(g-s*z)*(g-s*w)/(g*(g-d)); kd=(d-s*z)*(d-s*w)/(d*(d-g))
        add((b+1,i,j), Wt*s**2*t**(-i-j)*(1-t**(b+1)))
        for n in range(b+1):
            add((n,i+1,j), -Wt*kg*g**(n-b-1))
            add((n,i,j),   Wt*t**n*(kg*g**(n-b-1)+kd*d**(n-b-1)))
            add((n,i,j+1), -Wt*kd*d**(n-b-1))
    br=sum(t**l for l in range(k))
    return {key:sp.factor(sp.cancel(v/((1-t)*br))) for key,v in new.items() if sp.cancel(v)!=0}
def poch(x,N): return sp.prod([1-x*t**i for i in range(N)])
def al(J): return 0 if J<0 else sp.prod([s-t**i for i in range(1,J+1)])/poch(t,J)
def c(n,j): return 0 if (j<0 or j>n) else poch(s,n-j)/poch(t,n-j)*(al(j)-s*t**(n-j)*al(j-1))
def N(n,j): return t**(-n*j)*c(n,j)
def K(i,j): return sp.prod([(t**p*z-s*w)/(t**p*z-t**j*w) for p in range(i)])*sp.prod([(s*z-t**r*w)/(t**i*z-t**r*w) for r in range(j)])
def TC(k):
    out={}
    for b0 in range(k+1):
        for n1 in range(k-b0+1):
            n2=k-b0-n1
            for i in range(n1+1):
                for j in range(n2+1):
                    v=s**(2*b0+n1-i+n2-j)*t**(-b0*(i+j)-(n1-i)*j-(n2-j)*i)*K(i,j)*N(n1,i)*N(n2,j)*z**(-n1)*w**(-n2)
                    out[(b0,i,j)]=out.get((b0,i,j),0)+v
    return out
old=pickle.load(open('../day208/gf_Gamma_upto3.pkl','rb'))
K_MAX=int(sys.argv[1]) if len(sys.argv)>1 else 4
G={(0,0,0):sp.Integer(1)}; ok=True
for k in range(1,K_MAX+1):
    G=step(G,k)
    tc=TC(k)
    keys=set(G)|set(tc)
    d1=all(sp.cancel(G.get(q,0)-tc.get(q,0))==0 for q in keys)
    msg=f'k={k}: step==TC {d1}'
    if k in old:
        o={ (lam[0] if lam else 0,i,j):v for (lam,i,j),v in old[k].items()}
        assert all(len(lam)<=1 for (lam,i,j) in old[k])
        keys2=set(G)|set(o)
        d2=all(sp.cancel(G.get(q,0)-o.get(q,0))==0 for q in keys2)
        msg+=f'; step==day208 gf_recursion {d2}'; ok&=d2
    ok&=d1; print(msg,flush=True)
print('ALL OK' if ok else 'FAIL')

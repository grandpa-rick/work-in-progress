"""Day 212: closed ell-column step map vs (a) day210 gf_recursion_ell pickles, (b) (★ℓ) formula T_k.
Step map on W e_b prod_c E(gam_c), gam_c=t^{i_c} z_c, output x(1-t):
  e_{b+1}: W c (1-t^{b+1}),  c = s^ell t^{-S}
  n<=b: e_n: W[-sum_c kap_c gam_c^{n-b-1}(E(t gam_c) - t^n E(gam_c))]
  kap_c = prod_{c'}(gam_c - s z_c')/(gam_c prod_{c'!=c}(gam_c-gam_c'))."""
import sympy as sp, pickle, sys, itertools
s,t=sp.symbols('s t')
def step(G,k,Z):
    ell=len(Z); new={}
    def add(key,v): new[key]=new.get(key,0)+v
    for (b,I),Wt in G.items():
        gam=[t**I[c]*Z[c] for c in range(ell)]
        S=sum(I); cc=s**ell*t**(-S)
        kap=[sp.prod([gam[c]-s*Z[d] for d in range(ell)])/(gam[c]*sp.prod([gam[c]-gam[d] for d in range(ell) if d!=c])) for c in range(ell)]
        add((b+1,I),Wt*cc*(1-t**(b+1)))
        for n in range(b+1):
            for c in range(ell):
                J=list(I); J[c]+=1; J=tuple(J)
                add((n,J),-Wt*kap[c]*gam[c]**(n-b-1))
                add((n,I), Wt*t**n*kap[c]*gam[c]**(n-b-1))
    br=sum(t**l for l in range(k))
    out={}
    for key,v in new.items():
        v=sp.cancel(v/((1-t)*br))
        if v!=0: out[key]=v
    return out
def poch(x,N): return sp.prod([1-x*t**i for i in range(N)])
def al(J): return 0 if J<0 else sp.prod([s-t**i for i in range(1,J+1)])/poch(t,J)
def cf(n,j): return 0 if (j<0 or j>n) else poch(s,n-j)/poch(t,n-j)*(al(j)-s*t**(n-j)*al(j-1))
def N(n,j): return t**(-n*j)*cf(n,j)
def Kp(i,j,z,w): return sp.prod([(t**p*z-s*w)/(t**p*z-t**j*w) for p in range(i)])*sp.prod([(s*z-t**r*w)/(t**i*z-t**r*w) for r in range(j)])
def comps(n,ell):
    if ell==1: yield (n,); return
    for a in range(n+1):
        for r in comps(n-a,ell-1): yield (a,)+r
def V(n,I,Z):
    ell=len(Z); tot=0
    KI=sp.prod([Kp(I[a],I[b],Z[a],Z[b]) for a in range(ell) for b in range(a+1,ell)])
    for ns in comps(n,ell):
        if any(ns[c]<I[c] for c in range(ell)): continue
        tot+=s**((ell-1)*sum(ns[c]-I[c] for c in range(ell)))*t**(-sum((ns[c]-I[c])*I[d] for c in range(ell) for d in range(ell) if d!=c)) \
            *sp.prod([N(ns[c],I[c])*Z[c]**(-ns[c]) for c in range(ell)])
    return KI*tot
def T(k,Z):
    ell=len(Z); out={}
    for b in range(k+1):
        for I in itertools.product(range(k-b+1),repeat=ell):
            if sum(I)>k-b: continue
            v=(s**ell*t**(-sum(I)))**b*V(k-b,I,Z)
            if sp.cancel(v)!=0: out[(b,I)]=v
    return out
if __name__=='__main__':
    ell=int(sys.argv[1]); K=int(sys.argv[2]); Z=sp.symbols(f'z1:{ell+1}')
    old=None
    try:
        old=pickle.load(open(f'../day210/gamma_ell{ell}_k{sys.argv[3]}.pkl','rb')) if len(sys.argv)>3 else None
    except Exception as ex: print('no pickle',ex)
    G={(0,(0,)*ell):sp.Integer(1)}; ok=True
    for k in range(1,K+1):
        G=step(G,k,Z); tk=T(k,Z)
        d1=all(sp.cancel(G.get(q,0)-tk.get(q,0))==0 for q in set(G)|set(tk))
        msg=f'ell={ell} k={k}: step==T {d1} ({len(G)} terms)'; ok&=d1
        if old is not None and k in old:
            oz=old[k]; Zo=None
            o={}
            for (lam,I),v in oz.items():
                assert len(lam)<=1; o[(lam[0] if lam else 0,I)]=v
            # pickled symbols z1..zell same names
            d2=all(sp.cancel(G.get(q,0)-o.get(q,0))==0 for q in set(G)|set(o))
            msg+=f'; step==day210 gf_recursion_ell {d2}'; ok&=d2
        print(msg,flush=True)
    print('ALL OK' if ok else 'FAIL')

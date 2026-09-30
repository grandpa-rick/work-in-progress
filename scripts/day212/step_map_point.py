"""Day 212: step map iterated == T_k at exact rational points (s,t,z_c in Q). Larger (ell,k)."""
import sys, itertools
from fractions import Fraction as Fr
def run(ell,K,s,t,Z):
    def poch(x,N):
        r=Fr(1)
        for i in range(N): r*=1-x*t**i
        return r
    def al(J):
        if J<0: return Fr(0)
        r=Fr(1)
        for i in range(1,J+1): r*=s-t**i
        return r/poch(t,J)
    def cf(n,j): return Fr(0) if (j<0 or j>n) else poch(s,n-j)/poch(t,n-j)*(al(j)-s*t**(n-j)*al(j-1))
    def N(n,j): return t**(-n*j)*cf(n,j)
    def Kp(i,j,z,w):
        r=Fr(1)
        for p in range(i): r*=(t**p*z-s*w)/(t**p*z-t**j*w)
        for q in range(j): r*=(s*z-t**q*w)/(t**i*z-t**q*w)
        return r
    def comps(n,l):
        if l==1: yield (n,); return
        for a in range(n+1):
            for r in comps(n-a,l-1): yield (a,)+r
    def V(n,I):
        KI=Fr(1)
        for a in range(ell):
            for b in range(a+1,ell): KI*=Kp(I[a],I[b],Z[a],Z[b])
        tot=Fr(0)
        for ns in comps(n,ell):
            if any(ns[c]<I[c] for c in range(ell)): continue
            v=s**((ell-1)*sum(ns[c]-I[c] for c in range(ell)))*t**(-sum((ns[c]-I[c])*I[d] for c in range(ell) for d in range(ell) if d!=c))
            for c in range(ell): v*=N(ns[c],I[c])*Z[c]**(-ns[c])
            tot+=v
        return KI*tot
    def T(k):
        out={}
        for b in range(k+1):
            for I in itertools.product(range(k-b+1),repeat=ell):
                if sum(I)>k-b: continue
                v=(s**ell*t**(-sum(I)))**b*V(k-b,I)
                if v!=0: out[(b,I)]=v
        return out
    def step(G,k):
        new={}
        def add(key,v): new[key]=new.get(key,0)+v
        for (b,I),W in G.items():
            gam=[t**I[c]*Z[c] for c in range(ell)]; cc=s**ell*t**(-sum(I))
            kap=[]
            for c in range(ell):
                num=Fr(1); den=gam[c]
                for d in range(ell): num*=gam[c]-s*Z[d]
                for d in range(ell):
                    if d!=c: den*=gam[c]-gam[d]
                kap.append(num/den)
            add((b+1,I),W*cc*(1-t**(b+1)))
            for n in range(b+1):
                for c in range(ell):
                    J=list(I); J[c]+=1; J=tuple(J)
                    add((n,J),-W*kap[c]*gam[c]**(n-b-1)); add((n,I),W*t**n*kap[c]*gam[c]**(n-b-1))
        br=sum(t**l for l in range(k))
        return {q:v/((1-t)*br) for q,v in new.items() if v!=0}
    G={(0,(0,)*ell):Fr(1)}; ok=True
    for k in range(1,K+1):
        G=step(G,k); tk=T(k)
        d=all(G.get(q,0)==tk.get(q,0) for q in set(G)|set(tk))
        print(f'ell={ell} k={k}: step==T {d} ({len(tk)} terms)',flush=True); ok&=d
    return ok
ok=True
pts=[(Fr(3,7),Fr(5,11),[Fr(2,3),Fr(9,8),Fr(13,5),Fr(17,19),Fr(4,29)]),(Fr(-2,9),Fr(7,4),[Fr(5,2),Fr(-3,7),Fr(11,6),Fr(1,13),Fr(8,3)])]
for (s,t,Z) in pts:
    for ell,K in [(1,6),(2,6),(3,6),(4,5),(5,4)]:
        ok&=run(ell,K,s,t,Z[:ell])
print('ALL OK' if ok else 'FAIL')

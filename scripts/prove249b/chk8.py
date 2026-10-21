# sanity: with no P-branches firing, iota options must reproduce psi_opts0 lim=1 options exactly
import scert8 as S, itertools
S.NV=9; bad=0; tot=0
for r in range(1,4):
    for lam in itertools.combinations(range(1,8),r):
        for th in itertools.product(range(-3,4),repeat=r):
            S.LEMMAP=False; a=sorted(S.psi_opts0(lam,th,8,1),key=str)
            S.LEMMAP=True; S.IOTA=True
            # compare only when the lim0 call on thv uses the generic branch: emulate by disabling P inside
            thv=tuple(lam[j]-2*(j+1)+1-th[j] for j in range(r))
            S.LEMMAP=False
            S.LEMMAP=True
            saveP=(S.LEMMAP1,S.LEMMAP2)
            import types
            # temporarily force generic lim0 branch
            orig0=S.psi_opts0
            def gen0(l,t,n,lim,orig0=orig0):
                S.LEMMAP=False; o=orig0(l,t,n,lim); S.LEMMAP=True; return o
            S.psi_opts0=lambda l,t,n,lim: gen0(l,t,n,lim) if lim==0 else orig0(l,t,n,lim)
            b=sorted(S.psi_opts(lam,th,8,1),key=str)
            S.psi_opts0=orig0
            tot+=1
            norm=lambda L:sorted([(o,s,tuple(sorted((k,v) for k,v in d.items() if v)),w,nl,nt) for o,s,d,w,nl,nt in L])
            if norm(a)!=norm(b): bad+=1; 
            if norm(a)!=norm(b) and bad<3: print(lam,th,norm(a),norm(b))
print('generic iota vs direct lim=1:',tot,'cases, mismatches',bad)

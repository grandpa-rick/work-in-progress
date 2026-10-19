from mp import *
for c in (-3,2,5):
  for n,a in [(4,1),(4,2),(5,2),(5,3)]:
    sA=Sys(n,11,a,n); sB=Sys(n-1,11)
    rat=set()
    for l,m,nu in [((),(),tuple(range(1,n+1))),((),(1,),tuple(range(1,n+1))),((1,),(2,),(1,2,3,n)),((),(2,),(2,3,n)),((1,),(1,),(a,3,n) if a<3 else (1,a,n)),((),(),(a,n))]:
        if not EA(n,a,n,l,m,nu): continue
        r1=A(sA,a,n,l,m,nu,c); r2=B(sB,a,l,m,nu[:-1],c)
        rat.add(r1*inv(r2)%P if r2 else ('B0',r1))
    sB2=Sys(n,11); sC=Sys(n-1,11); rat2=set()
    for l,m,k in [((),(),tuple(range(1,n+1))),((),(1,),(1,3,n)),((1,),(2,),(2,3,n)),((2,),(2,),(1,2,n)),((),(),(n,))]:
        if not EB(n,n,l,m,k): continue
        r1=B(sB2,n,l,m,k,c); r2=C(sC,l,m,k[:-1],c-2)
        rat2.add(r1*inv(r2)%P if r2 else ('C0',r1))
    print(c,n,a,'A->B ratios',len(rat),[r for r in rat if isinstance(r,tuple)],' B->C(c-2) ratios',len(rat2),[r for r in rat2 if isinstance(r,tuple)])

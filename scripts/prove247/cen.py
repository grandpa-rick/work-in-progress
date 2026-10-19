import sys,time; from mp import *
typ=sys.argv[1]; n=int(sys.argv[2]); cs=list(map(int,sys.argv[3].split(',')))
Ss=subsets(n); t0=time.time()
sy=[Sys(n,sd) for sd in (1,2)]
for c in cs:
    zs=[];tot=0
    for i,l in enumerate(Ss):
        for m in Ss[i:]:
            for k in Ss:
                if typ=='C':
                    if not EC(n,l,m,k): continue
                    tot+=1
                    if all(C(s,l,m,k,c)==0 for s in sy): zs.append((l,m,k))
                elif typ=='B':
                    for a in range(1,n+1):
                        if not EB(n,a,l,m,k): continue
                        tot+=1
                        if all(B(s,a,l,m,k,c)==0 for s in sy): zs.append((a,l,m,k))
    print(typ,n,'c',c,'tot',tot,'zeros',len(zs),zs[:8],'%.0fs'%(time.time()-t0),flush=True)

import sys,time,itertools; from mp import *
n=int(sys.argv[1]); cs=list(map(int,sys.argv[2].split(','))); full=len(sys.argv)>3
Ss=subsets(n); t0=time.time()
pairs=list(itertools.combinations(range(1,n+1),2))
sy={p:[Sys(n,sd,p[0],p[1]) for sd in (1,2)] for p in pairs}
for c in cs:
    zs=[];tot=0
    for i,l in enumerate(Ss):
        for m in Ss[i:]:
            for k in Ss:
                for a,b in pairs:
                    if not EA(n,a,b,l,m,k): continue
                    tot+=1
                    if all(A(s,a,b,l,m,k,c)==0 for s in sy[a,b]): zs.append((a,b,l,m,k))
    print('A',n,'c',c,'tot',tot,'zeros',len(zs),zs if full else zs[:8],'%.0fs'%(time.time()-t0),flush=True)

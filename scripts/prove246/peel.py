import sys, itertools
sys.path.insert(0,'.')
from fast import subsets, geq
def E(l,m,nu,a,b,n):
    rest=[c for c in range(1,n+1) if c not in (a,b)]; out=[]
    for k in range(len(rest)+1):
        for Mpp in itertools.combinations(rest,k):
            M=tuple(sorted(Mpp+(a,b)))
            if geq(nu,M) and geq(Mpp,l) and geq(Mpp,m): out.append(Mpp)
    return out
def minus(k): return k[:-1]
def peel(l,m,nu,a,b):
    steps=0
    while True:
        n=max(nu); Es=E(l,m,nu,a,b,n)
        assert Es
        if n==b: return ('ii',l,m,nu,Es,steps)
        if not any(n in M for M in Es): return ('iii',l,m,nu,Es,steps)
        l,m,nu=minus(l),minus(m),minus(nu); steps+=1
        if len(nu)<2: return ('small',l,m,nu,E(l,m,nu,a,b,max(nu) if nu else 0) if nu else [],steps)
N=int(sys.argv[1]); Ss=subsets(N)
tri=[(l,m,nu) for i,l in enumerate(Ss) for m in Ss[i:] for nu in Ss if geq(nu,l) and geq(nu,m) and len(nu)-max(len(l),len(m))>=2]
stat={}; ex={}
for a,b in itertools.combinations(range(1,N+1),2):
    for l,m,nu in tri:
        if not E(l,m,nu,a,b,N): continue
        r=peel(l,m,nu,a,b)
        key=(r[0], 'single' if len(r[4])==1 else 'multi')
        stat[key]=stat.get(key,0)+1; ex.setdefault(key,[]).append(((a,b),l,m,nu,'->',r[1],r[2],r[3],r[4]))
print(stat)
for k,v in ex.items():
    print(k); [print('   ',x) for x in v[:8]]

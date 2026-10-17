import itertools,sys
from engine import subsets, geq
def cover(N):
    S=subsets(N); tot=0; ct=cs=0; unc=[]
    for i,l in enumerate(S):
        for m in S[i:]:
            for nu in S:
                if not(geq(nu,l) and geq(nu,m)): continue
                if len(nu)-max(len(l),len(m))<2: continue
                tot+=1
                top=geq(nu[:-2],l) and geq(nu[:-2],m)
                single=False
                for a,b in itertools.combinations(nu,2):
                    pp=tuple(x for x in nu if x not in (a,b))
                    if not(geq(pp,l) and geq(pp,m)): continue
                    D=[M for M in S if geq(nu,M) and a in M and b in M and geq(tuple(x for x in M if x not in(a,b)),l) and geq(tuple(x for x in M if x not in(a,b)),m)]
                    if D==[nu]: single=True;break
                ct+=top; cs+=single
                if not(top or single): unc.append((l,m,nu))
    return tot,ct,cs,unc
for N in range(2,int(sys.argv[1])+1):
    tot,ct,cs,unc=cover(N)
    print('N',N,'jump>=2 triples (lam<=mu order)',tot,'top-two',ct,'single-defect',cs,'uncovered',len(unc))
    if N==4: print(unc)

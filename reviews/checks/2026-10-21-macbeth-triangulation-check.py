"""Rick's checks for MacBeth 'weld-obstruction triangulation' (UID 425). Exact arithmetic where it matters."""
from fractions import Fraction as Fr
import itertools, random
import numpy as np
random.seed(20261021); rng=np.random.default_rng(20261021)

def K(M):  # row normalization, identity on zero rows
    M=np.array(M,dtype=object); out=M.copy()
    for i in range(M.shape[0]):
        r=sum(M[i]); 
        if r!=0: out[i]=[x/r for x in M[i]]
    return out
def Kcol(M): return K(np.array(M,dtype=object).T).T   # breaking control
def mm(A,B): return np.array(A,dtype=object).dot(np.array(B,dtype=object))
def eq(A,B): return all(x==y for x,y in zip(np.ravel(A),np.ravel(B)))

print("== C1 Lemma 9.3 witness (exact)")
f=[[Fr(2,3),Fr(1,3)]]; g=[[Fr(2,3),Fr(1,3),0],[0,Fr(1,2),Fr(1,2)]]; h=[[1,0],[0,0],[0,1]]
gh=mm(g,h); print(" gh =",gh.tolist()," masses",[sum(r) for r in gh])
L=K(mm(mm(f,g),h)); R=K(mm(f,K(gh))); LL=K(mm(K(mm(f,g)),h))
print(" (f#g)#h =",LL.tolist()," K(fgh) =",L.tolist()," f#(g#h) =",R.tolist())
print(" D_gh =",[1/sum(r) for r in gh])

def rand_sub(n,m,zero_p=0.25,exact=True,normalized=False):
    M=[]
    for i in range(n):
        if random.random()<zero_p: M.append([Fr(0)]*m); continue
        w=[Fr(random.randint(0,4)) for _ in range(m)]
        if sum(w)==0: w[random.randrange(m)]=Fr(1)
        s=sum(w); scale=Fr(1) if normalized else Fr(random.randint(1,5),5)
        M.append([x/s*scale for x in w])
    return M

print("== C2 (LA) for row-normalization, all substochastic x,c; control: column normalization")
bad=0; badc=0; N=3000
for t in range(N):
    a,b,c_=random.randint(1,4),random.randint(1,4),random.randint(1,4)
    x=rand_sub(a,b); c=rand_sub(b,c_)
    if not eq(K(mm(K(x),c)),K(mm(x,c))): bad+=1
    if not eq(Kcol(mm(Kcol(x),c)),Kcol(mm(x,c))): badc+=1
print(f" LA violations row-K {bad}/{N}; column-K (control) {badc}/{N}")

print("== C3 Thm 8.1 and defect census on Norm; supports of the two bracketings")
nd=0; supp_diff=0; comparable_when_diff=0; N=3000
for t in range(N):
    a,b,c_,d=[random.randint(1,4) for _ in range(4)]
    f=rand_sub(a,b,normalized=True); g=rand_sub(b,c_,normalized=True); h=rand_sub(c_,d,normalized=True)
    Lb=K(mm(K(mm(f,g)),h)); assert eq(Lb,K(mm(mm(f,g),h)))
    Rb=K(mm(f,K(mm(g,h))))
    if not eq(Lb,Rb):
        nd+=1
        le=all(x<=y for x,y in zip(np.ravel(Lb),np.ravel(Rb))); ge=all(x>=y for x,y in zip(np.ravel(Lb),np.ravel(Rb)))
        if le or ge: comparable_when_diff+=1
    if any((x>0)!=(y>0) for x,y in zip(np.ravel(Lb),np.ravel(Rb))): supp_diff+=1
print(f" Thm 8.1 OK on {N}; non-assoc triples {nd}; pointwise-comparable among those {comparable_when_diff}; support mismatches {supp_diff}/{N}")

print("== C4 recovery regime: h total (g arbitrary Norm, with bot rows) => associative; control h with a bot row")
na_total=0; na_ctrl=0; N=2000
for t in range(N):
    a,b,c_,d=[random.randint(1,4) for _ in range(4)]
    f=rand_sub(a,b,normalized=True); g=rand_sub(b,c_,zero_p=0.4,normalized=True); h=rand_sub(c_,d,zero_p=0,normalized=True)
    if not eq(K(mm(mm(f,g),h)),K(mm(f,K(mm(g,h))))): na_total+=1
    c3=max(c_,3); g2=rand_sub(b,c3,zero_p=0,normalized=True); h2=rand_sub(c3,max(d,2),zero_p=0,normalized=True); h2[0]=[Fr(0)]*len(h2[0])
    if not eq(K(mm(mm(f,g2),h2)),K(mm(f,K(mm(g2,h2))))): na_ctrl+=1
print(f" h total: non-assoc {na_total}/{N} (g has bot rows ~40%); control (h has bot row, mid>=3): non-assoc {na_ctrl}/{N}")

print("== C5 minimal shape a->b->c->d for a Norm defect (random search 600 per shape)")
for shape in [(1,2,2,2),(1,2,3,1),(1,1,3,2),(1,2,3,2),(2,2,2,3),(1,3,2,3),(1,2,4,2)]:
    a,b,c_,d=shape; cnt=0
    for t in range(600):
        f=rand_sub(a,b,zero_p=0.1,normalized=True); g=rand_sub(b,c_,zero_p=0.1,normalized=True); h=rand_sub(c_,d,zero_p=0.35,normalized=True)
        if not eq(K(mm(mm(f,g),h)),K(mm(f,K(mm(g,h))))): cnt+=1
    print(f"  shape {shape}: defects {cnt}/600")

print("== C6 monoids: normalizers with (K1),(U),(LA), # non-assoc")
def monoids(n):
    # element 0 = identity
    cells=[(i,j) for i in range(1,n) for j in range(1,n)]
    for vals in itertools.product(range(n),repeat=len(cells)):
        T=[[0]*n for _ in range(n)]
        for i in range(n): T[0][i]=i; T[i][0]=i
        for (i,j),v in zip(cells,vals): T[i][j]=v
        if all(T[T[x][y]][z]==T[x][T[y][z]] for x in range(n) for y in range(n) for z in range(n)): yield T
def handed(T,Kmap,useU=True):
    n=len(T); Km=lambda x:Kmap[x]
    if any(Km(Km(x))!=Km(x) for x in range(n)): return None
    if useU and Km(0)!=0: return None
    if any(Km(T[Km(x)][c])!=Km(T[x][c]) for x in range(n) for c in range(n)): return None
    Im=[x for x in range(n) if Km(x)==x]; sh=lambda a,b:Km(T[a][b])
    return any(sh(sh(a,b),c)!=sh(a,sh(b,c)) for a in Im for b in Im for c in Im)
for n in [1,2,3,4]:
    ms=list(monoids(n)); hits=0; closed_hits=0
    for T in ms:
        for Kmap in itertools.product(range(n),repeat=n):
            r=handed(T,Kmap)
            if r:
                hits+=1
                Im=[x for x in range(n) if Kmap[x]==x]
                if all(T[a][b] in Im for a in Im for b in Im): closed_hits+=1
    print(f"  order {n}: {len(ms)} monoid tables (labelled, e=0), handed normalizers {hits}, of which Im K closed under . : {closed_hits}")

print("== C7 the Lean model: M={e,z,a,b}, z zero, aa=z, ab=a, bb=b, K(z)=b; b.a undetermined in the PDF")
e,z,a,b=0,1,2,3
for ba,name in [(z,'z'),(a,'a')]:
    T=[[0]*4 for _ in range(4)]
    for i in range(4): T[e][i]=i; T[i][e]=i; T[z][i]=z if i!=e else z; T[i][z]=z
    T[e][z]=z; T[z][e]=z
    T[a][a]=z; T[a][b]=a; T[b][b]=b; T[b][a]=ba
    assoc=all(T[T[x][y]][w]==T[x][T[y][w]] for x in range(4) for y in range(4) for w in range(4))
    Kmap=[e,b,a,b]
    LA=all(Kmap[T[Kmap[x]][c]]==Kmap[T[x][c]] for x in range(4) for c in range(4))
    sh=lambda p,q:Kmap[T[p][q]]
    print(f"  b.a={name}: monoid {assoc}; LA {LA}; (a#a)#a={'ezab'[sh(sh(a,a),a)]}, a#(a#a)={'ezab'[sh(a,sh(a,a))]}")

print("== C8 principal normalizers K=k.x, x.k, k.x.k (k idempotent), WITHOUT (U): # associative? (orders<=4)")
viol=0; tot=0; withU_nontriv=0
for n in [2,3,4]:
    for T in monoids(n):
        for k in range(n):
            if T[k][k]!=k: continue
            for form in ['kx','xk','kxk']:
                Kmap=[T[k][x] if form=='kx' else T[x][k] if form=='xk' else T[T[k][x]][k] for x in range(n)]
                if any(Kmap[Kmap[x]]!=Kmap[x] for x in range(n)): continue
                Im=[x for x in range(n) if Kmap[x]==x]; sh=lambda p,q:Kmap[T[p][q]]
                tot+=1
                if any(sh(sh(p,q),r)!=sh(p,sh(q,r)) for p in Im for q in Im for r in Im): viol+=1
                if Kmap[0]==0 and any(Kmap[x]!=x for x in range(n)): withU_nontriv+=1
print(f"  principal normalizers tested {tot}; non-associative {viol}; satisfying (U) but K!=id {withU_nontriv}")

print("== C9 octonions as twisted group algebra of Z2^3 (Albuquerque-Majid): associator = coboundary of a 2-cochain")
G=list(itertools.product([0,1],repeat=3))
add=lambda x,y:tuple((p+q)%2 for p,q in zip(x,y))
def F(x,y):
    s=sum(x[i]*y[j] for i in range(3) for j in range(i,3)) + y[0]*x[1]*x[2]+x[0]*y[1]*x[2]+x[0]*x[1]*y[2]
    return (-1)**(s%2)
phi=lambda x,y,w: F(x,y)*F(add(x,y),w)*F(y,w)*F(x,add(y,w))  # +-1 so inverse = itself
nonassoc=sum(1 for x in G for y in G for w in G if phi(x,y,w)==-1)
cocyc=all(phi(y,w,u)*phi(x,add(y,w),u)*phi(x,y,w)==phi(add(x,y),w,u)*phi(x,y,add(w,u)) for x in G for y in G for w in G for u in G)
def omul(A,B):
    C={g:0 for g in G}
    for x in G:
        for y in G: C[add(x,y)]+=A[x]*B[y]*F(x,y)
    return C
def rv(): return {g:random.randint(-3,3) for g in G}
alt=all((lambda A,B:omul(omul(A,A),B)==omul(A,omul(A,B)))(rv(),rv()) for _ in range(200))
nm=lambda A:sum(v*v for v in A.values())
normmult=all((lambda A,B:nm(omul(A,B))==nm(A)*nm(B))(rv(),rv()) for _ in range(200))
print(f"  non-associative basis triples {nonassoc}/512; phi=dF is a 3-cocycle: {cocyc}; left-alternative: {alt}; norm multiplicative: {normmult}")

print("== C10 order-4 handed normalizers up to isomorphism (permutations fixing e)")
seen=set()
for T in monoids(4):
    for Kmap in itertools.product(range(4),repeat=4):
        if handed(T,Kmap):
            reps=[]
            for p in itertools.permutations([1,2,3]):
                P=[0]+list(p); Pi=[P.index(i) for i in range(4)]
                T2=tuple(tuple(P[T[Pi[i]][Pi[j]]] for j in range(4)) for i in range(4)); K2=tuple(P[Kmap[Pi[i]]] for i in range(4))
                reps.append((T2,K2))
            seen.add(min(reps))
print(f"  iso classes: {len(seen)}")
for T2,K2 in seen: print("  table",T2,"K",K2)

print("== C11 group transversals (Baer/Lal right loops) are normalizers in the sense of Def 7.1")
from itertools import permutations
G=list(permutations(range(3))); comp=lambda p,q: tuple(q[p[i]] for i in range(3))  # diagrammatic: p then q
e3=(0,1,2)
for Hgens in [[(1,0,2)],[(1,2,0)]]:
    H={e3}; 
    while True:
        new={comp(x,y) for x in H|set(Hgens) for y in H|set(Hgens)}|H
        if new==H: break
        H=new
    cosets={}
    for g in G: cosets.setdefault(frozenset(comp(h,g) for h in H),[]).append(g)   # H acting on the LEFT: Hg
    cl=list(cosets.values())
    cnt_tr=0; cnt_handed=0
    for choice in itertools.product(*cl):
        S=set(choice)
        if e3 not in S: continue
        cnt_tr+=1
        rep={g:[s for s in S if s in [c for c in cl if g in c][0]][0] for g in G}
        Km=lambda g:rep[g]
        assert all(Km(Km(g))==Km(g) for g in G) and Km(e3)==e3
        assert all(Km(comp(Km(x),c))==Km(comp(x,c)) for x in G for c in G)   # (LA)
        sh=lambda a,b:Km(comp(a,b))
        if any(sh(sh(a,b),c)!=sh(a,sh(b,c)) for a in S for b in S for c in S): cnt_handed+=1
    print(f"  |H|={len(H)}: normalized transversals {cnt_tr}, all satisfy (K1),(U),(LA); non-associative # : {cnt_handed}")

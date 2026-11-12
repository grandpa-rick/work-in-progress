"""Rick referee check of Clio Q584 (2026-11-11), written from her printed definitions only.
Lattice: point (a,b)=a e0+b e1; D[d] unit steps; class = d mod 3. Wall W={a=0}, H={a>=0}.
Triangles U(a,b)=(a,b),(a+1,b),(a,b+1); D(a,b)=(a+1,b),(a+1,b+1),(a,b+1). In H iff a>=0.
Part 1: enumerate connected regions (polyiamonds) in H, size<=NMAX, check Lemma 3, Lemma 4,
        Prop 5 (two bulk legs), Prop 6 (4 edges, no wall legs => rhombus).
Part 2: brute-force puzzle realisations of both sides of the 2-strand RE (Def 7) in a window,
        with K = any enumerated two-bulk-leg region and crossings = any 4-edge no-wall region.
        Positive controls: one-crossing piece; lesion allowing strands to turn in a crossing."""
import itertools, sys
D=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def verts(t):
    k,a,b=t
    return [(a,b),(a+1,b),(a,b+1)] if k=='U' else [(a+1,b),(a+1,b+1),(a,b+1)]
def edges(t):
    v=verts(t); return [frozenset((v[i],v[(i+1)%3])) for i in range(3)]
def eclass(e):
    p,q=tuple(e); d=(q[0]-p[0],q[1]-p[1]); return D.index(d)%3 if d in D else D.index((-d[0],-d[1]))%3
def nbrs(t):
    k,a,b=t
    return [('D',a,b),('D',a-1,b),('D',a,b-1)] if k=='U' else [('U',a,b),('U',a+1,b),('U',a,b+1)]
def bnd(R):
    cnt={}
    for t in R:
        for e in edges(t): cnt[e]=cnt.get(e,0)+1
    return [e for e,c in cnt.items() if c==1]
def is_wall(e): return all(p[0]==0 for p in e)
def directed_counts(R):
    # orient each boundary edge CCW w.r.t. its owning triangle (U is CCW as listed; D listed CCW too)
    n=[0,0,0]; B=set(bnd(R))
    for t in R:
        v=verts(t)
        for i in range(3):
            p,q=v[i],v[(i+1)%3]
            if frozenset((p,q)) in B:
                d=D.index((q[0]-p[0],q[1]-p[1])); n[d%3]+= 1 if d<3 else -1
    return n
def canon(R):
    mb=min(t[2] for t in R); return frozenset((k,a,b-mb) for k,a,b in R)
def enum(NMAX):
    seeds=[frozenset([('U',a,0)]) for a in range(NMAX+1)]+[frozenset([('D',a,0)]) for a in range(NMAX+1)]
    level={canon(s) for s in seeds}; allr=set(level)
    for n in range(2,NMAX+1):
        new=set()
        for R in level:
            for t in R:
                for u in nbrs(t):
                    if u[1]>=0 and u not in R: new.add(canon(R|{u}))
        level=new; allr|=new
    return allr
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 8
regs=enum(NMAX)
print(f'Part 1: {len(regs)} connected regions in H (size<={NMAX}, a-offset<= size, up to b-translation)')
viol3=mono=0; two=set(); four_bad=0; four=0; three_bulk_trap=None
def is_rhombus(R):
    if len(R)!=2: return False
    t1,t2=tuple(R); return t2 in nbrs(t1)
for R in regs:
    B=bnd(R); wl=[e for e in B if is_wall(e)]; bl=[e for e in B if not is_wall(e)]
    n=directed_counts(R)
    if not (n[0]==n[2]==-n[1]): viol3+=1
    if len({eclass(e) for e in bl})<=1: mono+=1
    if len(bl)==2: two.add((R,tuple(sorted(eclass(e) for e in bl)),len(wl)))
    if len(wl)==0 and len(B)==4:
        four+=1
        if not is_rhombus(R): four_bad+=1
print(' Lemma 3 closure violations:',viol3)
print(' regions with <=1 bulk class (Lemma 4):',mono)
print(' two-bulk-leg regions:',[(sorted(R),c,w) for R,c,w in two])
print(f' four-edge no-wall regions: {four}, non-rhombi among them: {four_bad}')
trap=frozenset([('U',0,0),('D',0,0),('U',0,1)])
B=bnd(trap); print(' trapezoid bulk-leg classes:',sorted(eclass(e) for e in B if not is_wall(e)),' wall legs:',sum(map(is_wall,B)))

# ---------- Part 2: realisation search ----------
K_shapes=[R for R,c,w in two]          # translate along wall
X_shapes=[R for R in regs if len(bnd(R))==4 and not any(map(is_wall,bnd(R)))]
WIN_A=range(0,5); WIN_B=range(-5,6)
def shift(R,da,db): return frozenset((k,a+da,b+db) for k,a,b in R)
Ks=set(); Xs=set()
for R in K_shapes:
    for db in WIN_B: Ks.add(shift(R,0,db))
for R in X_shapes:
    ma=min(t[1] for t in R)
    for da in WIN_A:
        for db in WIN_B:
            S=shift(R,da-ma,db)
            if all(t[1]>=0 for t in S) and not any(map(is_wall,bnd(S))): Xs.add(S)
Ks=[(R,frozenset(e for e in bnd(R) if not is_wall(e))) for R in Ks]; Xs=list(Xs)
if len(sys.argv)>2 and sys.argv[2]=='lesionK':
    # LESION: K may be ANY unit triangle in the window (bulk allowed), any 2 of its 3 edges as legs,
    # the third edge an auxiliary 'boundary' leg (trivalent vertex with one leg absorbed).
    Ks=[]
    for k in 'UD':
        for a in WIN_A:
            for b in WIN_B:
                t=(k,a,b); E=edges(t)
                for i in range(3): Ks.append((frozenset([t]),frozenset(E[j] for j in range(3) if j!=i)))
    print(' LESION K active:',len(Ks),'K choices')
print(f'\nPart 2: {len(Ks)} K placements, {len(Xs)} crossing placements')
def shared(R,S):
    s=set(bnd(R))&set(bnd(S)); return list(s)
def opposite(X,e,f):
    # in a 4-edge region: e,f opposite iff they share no vertex
    return e!=f and not (e&f)
def search(side, lesion=False):
    """side 'L': A: open-X1-X2-KA-open ; B: open-X1-KB-X2-open
       side 'R': A: open-KA-X1-X2-open ; B: open-X1-KB-X2-open"""
    sols=0; classes=set()
    for KA,LA in Ks:
      for KB,LB in Ks:
        if KA&KB: continue
        for X1 in Xs:
          if X1&KA or X1&KB: continue
          for X2 in Xs:
            if X2&X1 or X2&KA or X2&KB: continue
            e12=shared(X1,X2); 
            if len(e12)!=1: continue
            e1B=shared(X1,KB); eB2=shared(KB,X2)
            if len(e1B)!=1 or len(eB2)!=1: continue
            if side=='L':
                e2A=shared(X2,KA)
                if len(e2A)!=1 or shared(X1,KA): continue
                # X1: A-ends = (openA, e12); B-ends = (openB, e1B)
                b1=bnd(X1); oA=[e for e in b1 if e not in (e12[0],e1B[0])]
                if len(oA)!=2: continue
                # need: A ends opposite, B ends opposite => e12 opposite some open, e1B opposite other
                okX1 = lesion or any(opposite(X1,e12[0],o) and opposite(X1,e1B[0],[p for p in oA if p!=o][0]) for o in oA)
                b2=set(bnd(X2)); rest=list(b2-{e12[0],e2A[0],eB2[0]})
                okX2 = len(rest)==1 and (lesion or (opposite(X2,e12[0],e2A[0]) and opposite(X2,eB2[0],rest[0])))
                kb=LA
                okKA = e2A[0] in kb
                kbB=LB
                okKB = set(kbB)=={e1B[0],eB2[0]}
            else:
                eA1=shared(KA,X1)
                if len(eA1)!=1 or shared(X2,KA): continue
                b1=set(bnd(X1)); rest=list(b1-{eA1[0],e12[0],e1B[0]})
                okX1 = len(rest)==1 and (lesion or (opposite(X1,eA1[0],e12[0]) and opposite(X1,e1B[0],rest[0])))
                b2=bnd(X2); oo=[e for e in b2 if e not in (e12[0],eB2[0])]
                okX2 = len(oo)==2 and (lesion or any(opposite(X2,e12[0],o) and opposite(X2,eB2[0],[p for p in oo if p!=o][0]) for o in oo))
                kb=LA
                okKA = eA1[0] in kb
                kbB=LB
                okKB = set(kbB)=={e1B[0],eB2[0]}
            if okX1 and okX2 and okKA and okKB:
                sols+=1
    return sols
def one_crossing():
    """K_B, then B_hi x A_lo at X, then K_A:  B: open-KB-X-open ; A: open-X-KA-open"""
    n=0
    for KA,LA in Ks:
      for KB,LB in Ks:
        if KA&KB: continue
        for X in Xs:
          if X&KA or X&KB: continue
          eB=shared(KB,X); eA=shared(X,KA)
          if len(eB)!=1 or len(eA)!=1 or eB[0] not in LB or eA[0] not in LA: continue
          b=bnd(X); rest=[e for e in b if e not in (eB[0],eA[0])]
          if len(rest)!=2: continue
          if any(opposite(X,eB[0],r) and opposite(X,eA[0],[p for p in rest if p!=r][0]) for r in rest): n+=1
    return n
print(' LHS realisations:',search('L'))
print(' RHS realisations:',search('R'))
print(' control: one-crossing pictures (should be >0):',one_crossing())
print(' lesion (strands may turn inside crossings) LHS/RHS (should be >0 if search is live):',search('L',True),search('R',True))

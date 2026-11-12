"""Fast adjacency-driven version of the Def-7 realisation search (Rick, referee of Clio Q584).
Vertices: crossings = unit rhombi (Prop 6, re-verified in q584_check.py); K = region with exactly
the two designated legs. Modes:
  wall   : K = wall triangle U(0,b) with its two bulk legs (Prop 5, re-verified)      [Clio's setting]
  lesionK: K = ANY unit triangle in H, any 2 of its 3 edges as legs (third = absorbed aux leg)
  lesionT: wall K, but strands may turn inside a crossing (ends need not be opposite)
Arcs = directly shared edges (Clio Def 7, literal)."""
import sys, itertools
D=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def verts(t):
    k,a,b=t
    return [(a,b),(a+1,b),(a,b+1)] if k=='U' else [(a+1,b),(a+1,b+1),(a,b+1)]
def edges(t):
    v=verts(t); return [frozenset((v[i],v[(i+1)%3])) for i in range(3)]
def nbrs(t):
    k,a,b=t
    return [('D',a,b),('D',a-1,b),('D',a,b-1)] if k=='U' else [('U',a,b),('U',a+1,b),('U',a,b+1)]
def bnd(R):
    c={}
    for t in R:
        for e in edges(t): c[e]=c.get(e,0)+1
    return frozenset(e for e,n in c.items() if n==1)
def wall(e): return all(p[0]==0 for p in e)
def ecl(e):
    p,q=tuple(e); d=(q[0]-p[0],q[1]-p[1]); return (D.index(d) if d in D else D.index((-d[0],-d[1])))%3
mode=sys.argv[1] if len(sys.argv)>1 else 'wall'
A=range(0,6); B=range(-6,7)
tris=[(k,a,b) for k in 'UD' for a in A for b in B]
X=[]
for t in tris:
    for u in nbrs(t):
        if u[1]>=0 and t<u:
            R=frozenset([t,u]); bd=bnd(R)
            if not any(map(wall,bd)): X.append((R,bd))
K=[]
for t in tris:
    E=edges(t)
    if mode=='lesionK':
        for i in range(3): K.append((frozenset([t]),frozenset(E[j] for j in range(3) if j!=i)))
    elif t[0]=='U' and t[1]==0:
        K.append((frozenset([t]),frozenset(e for e in E if not wall(e))))
turn = mode=='lesionT'
def opp(e,f): return e!=f and not (e&f)
def okX(bd,pairA,pairB):
    # pairA/pairB: the two arc-ends of each strand at this crossing (None = open end)
    ends=list(pairA)+list(pairB); named=[e for e in ends if e is not None]
    if len(set(named))!=len(named) or not set(named)<=bd: return False
    free=list(bd-set(named)); nfree=sum(e is None for e in ends)
    if len(free)!=nfree: return False
    for perm in itertools.permutations(free):
        it=iter(perm); fill=[e if e is not None else next(it) for e in ends]
        if turn or (opp(fill[0],fill[1]) and opp(fill[2],fill[3])): return True
    return False
Kby={}
for R,L in K:
    for e in L: Kby.setdefault(e,[]).append((R,L))
def count(side):
    n=0; cls=set()
    for (X1,b1),(X2,b2) in itertools.permutations(X,2):
        if X1&X2: continue
        s=b1&b2
        if len(s)!=1: continue
        e12=next(iter(s))
        for e1B in b1-{e12}:
          for eB2 in b2-{e12}:
            for KB,LB in Kby.get(e1B,[]):
              if LB!={e1B,eB2} or KB&X1 or KB&X2: continue
              if side=='L':   # A: open-X1-X2-KA-open ; B: open-X1-KB-X2-open
                for e2A in b2-{e12,eB2}:
                  for KA,LA in Kby.get(e2A,[]):
                    if KA&(X1|X2|KB): continue
                    if not okX(b1,(None,e12),(None,e1B)): continue
                    if not okX(b2,(e12,e2A),(eB2,None)): continue
                    n+=1; cls.add(tuple(ecl(e) for e in (e12,e1B,eB2,e2A)))
              else:           # A: open-KA-X1-X2-open ; B: open-X1-KB-X2-open
                for eA1 in b1-{e12,e1B}:
                  for KA,LA in Kby.get(eA1,[]):
                    if KA&(X1|X2|KB): continue
                    if not okX(b1,(eA1,e12),(None,e1B)): continue
                    if not okX(b2,(e12,None),(eB2,None)): continue
                    n+=1; cls.add(tuple(ecl(e) for e in (e12,e1B,eB2,eA1)))
    return n,sorted(cls)
print('mode',mode,'| #X',len(X),'#K',len(K))
for s in 'LR':
    n,c=count(s); print(f' side {s}: realisations {n}; class patterns (A-arc X1X2, B X1-KB, B KB-X2, A K-arc): {c[:6]}')

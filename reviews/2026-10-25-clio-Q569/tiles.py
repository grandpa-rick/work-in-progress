"""Independent re-derivation (Rick, 2026-10-25) of Clio Q569 Thm 4, Cor 5, Prop 12, 13, 16, Thm 9
from the step sequences printed in her PDF section 2. Integer coords (a,b) <-> a e0 + b e1."""
from fractions import Fraction as F
D=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
def add(p,q): return (p[0]+q[0],p[1]+q[1])
TILES={
 'tri0U':[(0,0),(2,0),(4,0)],'tri0D':[(1,0),(3,0),(5,0)],
 'tri1U':[(0,1),(2,1),(4,1)],'tri1D':[(1,1),(3,1),(5,1)],
 'rh0':[(0,1),(1,0),(3,1),(4,0)],
 'K':[(1,1),(1,0),(3,1),(3,0),(5,1),(5,0)],
 'tri':[(0,1),(0,0),(2,1),(2,0),(4,1),(4,0)],
 'hex':[(1,0),(2,1),(3,0),(4,1),(5,0),(0,1)],
 'hexR':[(0,0),(1,1),(2,0),(3,1),(4,0),(5,1)],
 'trap':[(1,0),(2,1),(4,0),(4,1),(0,1)],
}
# rhombi "+2 rotations": rotate by 120, 240
def rot(seq,r): return [((d+r)%6,l) for d,l in seq]
TILES['rh1']=rot(TILES['rh0'],2); TILES['rh2']=rot(TILES['rh0'],4)
def walk(seq,start=(0,0)):
    pts=[start]; edges={}
    for d,l in seq:
        q=add(pts[-1],D[d]); edges[frozenset((pts[-1],q))]=l; pts.append(q)
    assert pts[-1]==start, 'not closed'
    return pts[:-1],edges
def euc(p): return (F(p[0])+F(p[1],2), p[1])  # y scaled by 2/sqrt3; affine, fine for incidence
def area2(pts):
    s=0
    for i in range(len(pts)):
        x1,y1=euc(pts[i]); x2,y2=euc(pts[(i+1)%len(pts)]); s+=x1*y2-x2*y1
    return s
def tri_pts(t):
    k,a,b=t
    return [(a,b),(a+1,b),(a,b+1)] if k=='U' else [(a+1,b),(a+1,b+1),(a,b+1)]
def centroid(pts):
    return (F(sum(p[0] for p in pts),len(pts)),F(sum(p[1] for p in pts),len(pts)))
def inside(pt,poly):  # convex or not: ray casting in sheared coords (affine ok)
    x,y=pt; c=False; n=len(poly)
    for i in range(n):
        (x1,y1),(x2,y2)=poly[i],poly[(i+1)%n]
        if (y1>y)!=(y2>y):
            xi=x1+(y-y1)*(x2-x1)/F(y2-y1)
            if xi>x: c=not c
    return c
def triangles(pts):
    A=min(p[0] for p in pts)-2;B=max(p[0] for p in pts)+2
    C=min(p[1] for p in pts)-2;E=max(p[1] for p in pts)+2
    out=[]
    for k in 'UD':
        for a in range(A,B):
            for b in range(C,E):
                if inside(centroid(tri_pts((k,a,b))),pts): out.append((k,a,b))
    return out
def labelled(seq,start=(0,0)):
    pts,ed=walk(seq,start); return frozenset(triangles(pts)),frozenset(ed.items())
def normalise(tile):
    """translate so lexicographically least triangle is at origin"""
    tris,ed=tile; k,a,b=min(tris,key=lambda t:(t[1],t[2],t[0]))
    sh=lambda p:(p[0]-a,p[1]-b)
    return (frozenset((t[0],t[1]-a,t[2]-b) for t in tris),
            frozenset((frozenset(sh(p) for p in e),l) for e,l in ed))
# group action on lattice points: rotation by 60 = linear map D_i->D_{i+1}; reflection r0: d->-d (axis dir 0)
def lin_rho(p):  # e0->e1, e1->D2
    a,b=p; return (a*0+b*(-1), a*1+b*1)
def lin_ref(p):  # D_d -> D_{-d}: e0->e0, e1->D5=(1,-1)
    a,b=p; return (a+b,-b)
def apply_seq(seq,g):
    # act on the walk by mapping its vertex points, rebuild region/edges
    pts,ed=walk(seq)
    gp=[g(p) for p in pts]
    if area2(gp)<0: gp=gp[::-1]
    ged={frozenset(g(p) for p in e):l for e,l in ed.items()}
    return frozenset(triangles(gp)),frozenset(ged.items())
def group():
    els=[]
    for r in range(6):
        def g(p,r=r):
            for _ in range(r): p=lin_rho(p)
            return p
        els.append(('rho%d'%r,g))
        def h(p,r=r):
            p=lin_ref(p)
            for _ in range(r): p=lin_rho(p)
            return p
        els.append(('ref%d'%r,h))
    return els
def swap(tile):
    tris,ed=tile; return tris,frozenset((e,1-l) for e,l in ed)
if __name__=='__main__':
    for nm,s in TILES.items():
        pts,_=walk(s); print(nm,'closes, 2*area(sheared)=',area2(pts),'#tri=',len(triangles(pts)))
    # Thm 4: hex region = vertex star of its centroid; labels keyed by orientation
    for nm in ['hex','hexR']:
        pts,ed=walk(TILES[nm]); c=centroid(pts); assert c[0].denominator==1 and c[1].denominator==1
        v=(int(c[0]),int(c[1])); assert sorted(pts)==sorted(add(v,d) for d in D)
        lab={}
        for e,l in ed.items():
            p,q=tuple(e); t=[v,p,q]
            # orientation: U if triangle is U(a,b) form
            ori=[k for k in 'UD' for a in range(v[0]-2,v[0]+2) for b in range(v[1]-2,v[1]+2) if sorted(tri_pts((k,a,b)))==sorted(t)]
            lab.setdefault(ori[0],set()).add(l)
        print(nm,'centre',v,'labels by orientation of owning star triangle:',lab)
    N={nm:normalise(labelled(s)) for nm,s in TILES.items()}
    KTW={N[x] for x in ['tri0U','tri0D','tri1U','tri1D','rh0','rh1','rh2']}
    print('\nCor5/Prop12/13: g, g(hex)==hex?, ==hexR?, g preserves KTW?, g∘swap preserves KTW?, g∘swap(hex)==hex?')
    for nm,g in group():
        gh=normalise(apply_seq(TILES['hex'],g))
        gK={normalise(apply_seq(TILES[x],g)) for x in ['tri0U','tri0D','tri1U','tri1D','rh0','rh1','rh2']}
        gKs={normalise(swap(t)) for t in gK}
        # orientation test
        up=apply_seq([(0,0),(2,0),(4,0)],g)[0]; ori=list(up)[0][0]
        print(nm,'up->',ori,gh==N['hex'],gh==N['hexR'],gK==KTW,gKs==KTW,normalise(swap(gh))==N['hex'])
    # Thm 9: trap = hex ∩ {a>=0} with v on wall
    v=(0,1); start=add(v,D[4])  # hex walk starts with dir1; find start so centre is v
    for st in [add(v,d) for d in D]:
        pts,ed=walk(TILES['hex'],st)
        if centroid(pts)==(v[0],v[1]): break
    ht=set(triangles(pts)); keep={t for t in ht if all(p[0]>=0 for p in tri_pts(t))}
    tp,te=walk(TILES['trap'],(1,0)); tt=set(triangles(tp))
    print('\nThm9: hex∩{a>=0} triangles == trap triangles:',keep==tt, sorted(tt))
    print(' retained hex edges with a>=0 (both ends):',{tuple(sorted(e)):l for e,l in ed.items() if all(p[0]>=0 for p in e)})
    print(' trap edges:',{tuple(sorted(e)):l for e,l in te.items()})
    # Prop 16: 1-edges per direction class
    for nm,s in TILES.items():
        c=[0,0,0]
        for d,l in s:
            if l>=1: c[d%3]+=1
        print('Prop16',nm,c)

# Re-check of MacBeth uid456 (mirror-LA v2, content 6df0043) from the printed statements.
from itertools import permutations, product
# (1) Thm 3.1 second separator as printed: M={e,a,b}, aa=b, ab=ba=a, bb=b; K(e)=e,K(a)=a,K(b)=a
E,A,B=0,1,2
T={(E,E):E,(E,A):A,(E,B):B,(A,E):A,(B,E):B,(A,A):B,(A,B):A,(B,A):A,(B,B):B}
m=lambda x,y:T[(x,y)]; K={E:E,A:A,B:A}
els=[E,A,B]
assoc=all(m(m(x,y),z)==m(x,m(y,z)) for x in els for y in els for z in els)
LA=all(K[m(K[x],c)]==K[m(x,c)] for x in els for c in els)
mLA=all(K[m(c,K[x])]==K[m(c,x)] for x in els for c in els)
img=[x for x in els if K[x]==x]; closed=all(m(f,g) in img for f in img for g in img)
nassoc=all(K[m(K[m(f,g)],h)]==K[m(f,K[m(g,h)])] for f in img for g in img for h in img)
print("sep2 printed: assoc",assoc,"LA",LA,"mLA",mLA,"closed",closed,"Norm assoc",nassoc)
# negative control: the old (pre-456) table bb=a
T2=dict(T); T2[(B,B)]=A
print("old sep2 table assoc:",all(T2[(T2[(x,y)],z)]==T2[(x,T2[(y,z)])] for x in els for y in els for z in els))
# (2) Rem 5.2 / 4.3 over S3xC2 and S4: theta_K == core-coset relation; Prop 5.1 (1)<=>(3); Rem 4.3 two forms
def compose(p,q): return tuple(q[p[i]] for i in range(len(p)))  # diagrammatic: p then q
def inv(p):
    r=[0]*len(p)
    for i,x in enumerate(p): r[x]=i
    return tuple(r)
def gen_group(gens):
    n=len(gens[0]); G={tuple(range(n))}; fr=list(G)
    while fr:
        nf=[]
        for g in fr:
            for s in gens:
                h=compose(g,s)
                if h not in G: G.add(h); nf.append(h)
        fr=nf
    return sorted(G)
def subgroups(G):
    out=set()
    for a in G:
        for b in G:
            out.add(frozenset(gen_group([a,b])))
    return out
def run(name,G):
    e=tuple(range(len(G[0])))
    stats=dict(pairs=0,p51=0,core=0,img=0,between=0,theta=0)
    for H in subgroups(G):
        if len(H)==len(G): continue
        Hs=set(H)
        cos={}; 
        for g in G: cos.setdefault(frozenset(compose(h,g) for h in H),[]).append(g)
        cl=[c for c in cos.values() if e not in c]
        core={x for x in H if all(compose(compose(inv(g),x),g) in Hs for g in G)}
        normal=len(core)==len(H)
        # theta_K vs core-coset relation (independent of T)
        for a in G:
            for b in G:
                th=all(compose(compose(x,compose(a,inv(b))),inv(x)) in Hs for x in G)
                if th != (compose(a,inv(b)) in core): stats['theta']+=1
        for choice in product(*cl):
            Tt=[e]+list(choice); Ts=set(Tt); stats['pairs']+=1
            rep={}
            for t in Tt:
                for h in H: rep[compose(h,t)]=t
            Kf=lambda x:rep[x]
            sharp=lambda f,g:Kf(compose(f,g))
            assoc=all(sharp(sharp(f,g),h)==sharp(f,sharp(g,h)) for f in Tt for g in Tt for h in Tt)
            # (3): K(gh) theta gh  <=> (K(gh))(gh)^-1 in core
            c3=all(compose(Kf(compose(g,h)),inv(compose(g,h))) in core for g in Tt for h in Tt)
            phi_core=all(compose(compose(g,h),inv(Kf(compose(g,h)))) in core for g in Tt for h in Tt)
            Tbar={frozenset(compose(c,t) for c in core) for t in Tt}
            imgsub=all(frozenset(compose(c,compose(s,t)) for c in core) in Tbar for s in Tt for t in Tt)
            closed=all(compose(g,h) in Ts for g in Tt for h in Tt)
            if assoc!=c3: stats['p51']+=1
            if assoc!=phi_core: stats['core']+=1
            if assoc!=imgsub: stats['img']+=1
            if assoc and not closed and not normal: stats['between']+=1
    print(name,stats)
S3xC2=gen_group([(1,0,2,3,4),(1,2,0,3,4),(0,1,2,4,3)])
S4=gen_group([(1,0,2,3),(1,2,3,0)])
run("S3xC2",S3xC2); run("S4",S4)

# generic: poset (elements 0..n-1, relation less(a,b)), D subset; F(N)=#order-preserving g:P->[0,N], g even on D
import sympy as sp, itertools, random, sys
x=sp.Symbol('x')
def ideals(n, less):
    below=[frozenset(b for b in range(n) if less(b,a)) for a in range(n)]
    res=[frozenset()]; seen={frozenset()}; fr=[frozenset()]
    while fr:
        new=[]
        for I in fr:
            for a in range(n):
                if a not in I and below[a]<=I:
                    J=I|{a}
                    if J not in seen: seen.add(J); new.append(J)
        res+=new; fr=new
    return res
def counter(n, less, D):
    J=ideals(n,less); tr=[(I,K,bool((K-I)&D)) for I in J for K in J if I<=K]
    full=frozenset(range(n))
    def F(N):
        vec={frozenset():1}
        for v in range(0,N+1):
            nv={}
            for I,K,hd in tr:
                c=vec.get(I)
                if not c or (hd and v%2): continue
                nv[K]=nv.get(K,0)+c
            vec=nv
        return vec.get(full,0)
    return F
def EO(n, less, D):
    F=counter(n,less,D); vals=[F(N) for N in range(2*n+4)]
    E=sp.expand(sp.interpolate([(N,vals[N]) for N in range(0,len(vals),2)],x))
    O=sp.expand(sp.interpolate([(N,vals[N]) for N in range(1,len(vals),2)],x))
    return E,O
def closure(n, rel):
    L=[[False]*n for _ in range(n)]
    for a,b in rel: L[a][b]=True
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if L[i][k] and L[k][j]: L[i][j]=True
    return lambda a,b: L[a][b]
if __name__=='__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
    for trial in range(40):
        n=random.randint(2,6)
        # element 0 = unique min; random extra relations respecting index order
        rel=[(0,i) for i in range(1,n)]+[(a,b) for a in range(1,n) for b in range(a+1,n) if random.random()<0.35]
        less=closure(n,rel)
        # D: chain containing 0
        others=[i for i in range(1,n)]; random.shuffle(others)
        D=[0]
        for i in sorted(others[:random.randint(0,n-1)]):
            if less(D[-1],i): D.append(i)
        E,O=EO(n,less,frozenset(D))
        e1=E.subs(x,-1)
        pred = None
        if len(D)>=2:
            d1=D[1]; pred=sp.Rational(1,2**sum(1 for a in range(n) if not (a==d1 or less(d1,a))))
        else: pred=sp.Rational(1,2**n)
        print(n, D, 'rel',[(a,b) for a in range(n) for b in range(n) if less(a,b)], 'E(-1)=',e1, 'pred',pred, 'OK' if e1==pred else 'XX')

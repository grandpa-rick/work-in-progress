from pq import *
m=int(sys.argv[1]); j=sp.Symbol('j')
for r in range(1, m+1):
    for k in itertools.combinations(range(1, m+1), r):
        S,E,O,p,q = polys(k)
        print(k, sp.factor(E.subs(x,2*j)))

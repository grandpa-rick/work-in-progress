from pq import *
m=int(sys.argv[1])
for r in range(1, m+1):
    for k in itertools.combinations(range(1, m+1), r):
        S,E,O,p,q = polys(k)
        sg=(-1)**sum(k); K=2**max(k)
        Eo=[sg*E.subs(x,-M)*K for M in range(1,16,2)]
        Fo=[Fint(S,M)*K for M in range(1,16,2)]
        print(k, 'odd M=1,3..15: 2^k1*E°=',Eo,' 2^k1*F°=',Fo)

from pq import *
m=int(sys.argv[1])
for r in range(1, m+1):
    for k in itertools.combinations(range(1, m+1), r):
        S,E,O,p,q = polys(k)
        Q = sp.expand(q*2**(max(k)+1))
        print(k, 'Q(-1..6)=', [Q.subs(x,v) for v in range(-1,7)], ' Q(-M) M=2..10', [Q.subs(x,-v) for v in range(2,11)])

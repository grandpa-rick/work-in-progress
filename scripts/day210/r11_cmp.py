import sympy as sp
from ds3 import elam
from tc_extract import s, t
def qi(n): return sum(t**i for i in range(n))
for r in range(1, 7):
    C = {(r,1,1): s**3, (r,2): s**2*(1-s)*qi(2), (r+1,1): s*(1-s)*(qi(r+1)+t*qi(r)+s), (r+2,): (1-s)**2*qi(r+2)*(qi(r+1)+s)}
    C2 = {}
    for mu,v in C.items():
        key = tuple(sorted(mu, reverse=True)); C2[key] = C2.get(key,0)+v
    E = elam((1, r, 1)) if r>=1 else None   # e_1*(e_r*e_1): TC k=1 outer
    bad = [m for m in set(C2)|set(E) if sp.cancel(C2.get(m,0)-E.get(m,0))!=0]
    print('r=',r,'TC-route e_1*(e_r*e_1) vs Day207 DS(r,1,1) closed form:', 'MATCH' if not bad else bad)

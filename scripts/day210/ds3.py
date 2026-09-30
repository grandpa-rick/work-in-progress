"""DS at length 3 as a corollary of 207b (e_k*e_r) + (TC) (e_k*(e_a e_b)).
e_lam^{(q,t)} = e_{l1}*e_{l2}*e_{l3} (commutative; Day 196 convention t^{-sum C(l_i,2)} prod e_{l_i}(Y).1).
Route: e_k*(e_x*e_y) with k chosen; inner via 207b, outer via TC/207b."""
import sympy as sp, sys, itertools
from tc_extract import coeffs, ek_er, s, t

def add(d, mu, v): d[mu] = d.get(mu, 0) + v
def star(k, F):
    out = {}
    for mu, c in F.items():
        if len(mu) == 1: G = ek_er(k, mu[0])
        elif len(mu) == 2: G = coeffs(k, mu[0], mu[1])
        else: raise ValueError
        for nu, d in G.items(): add(out, nu, c*d)
    return {m: sp.factor(v) for m, v in out.items() if sp.cancel(v) != 0}
def elam(order):
    k, x, y = order
    return star(k, ek_er(x, y))
def dominates(mu, lam):
    a = b = 0
    for i in range(max(len(mu), len(lam))):
        a += mu[i] if i < len(mu) else 0; b += lam[i] if i < len(lam) else 0
        if a < b: return False
    return True
def nstat(lam): return sum(i*l for i, l in enumerate(lam))
def parts3(n):
    return [(a,b,c) for a in range(n,0,-1) for b in range(a,0,-1) for c in range(b,0,-1) if a+b+c==n]

if __name__ == '__main__':
    NMAX = int(sys.argv[1])
    allok = True
    store = {}
    for n in range(3, NMAX+1):
        for lam in parts3(n):
            E = elam((lam[2], lam[0], lam[1]))      # e_{l3} * (e_{l1} * e_{l2})
            store[lam] = E
            supp_ok = all(dominates(mu, lam) for mu in E)
            lead_ok = sp.cancel(E.get(lam, 0) - s**nstat(lam)) == 0
            q1_ok = all(sp.cancel(v.subs(s, 1)) == 0 for mu, v in E.items() if mu != lam)
            full = all(mu in E for mu in [m for m in itertools.chain(*[[tuple(sorted(p, reverse=True))] for p in []])])
            nd = sum(1 for m in E); ndom = None
            ok = supp_ok and lead_ok and q1_ok
            allok &= ok
            print(lam, 'supp⊆up-set:', supp_ok, ' lead=s^n(λ):', lead_ok, ' offdiag|_{q=1}=0:', q1_ok, ' #supp', nd, flush=True)
            # second order (commutativity cross-check) when parts differ
            if len(set(lam)) > 1:
                E2 = elam((lam[0], lam[1], lam[2]))
                diff = [m for m in set(E)|set(E2) if sp.cancel(E.get(m,0)-E2.get(m,0)) != 0]
                print('    order (l1 outer) vs (l3 outer):', 'AGREE' if not diff else f'DIFFER {diff}', flush=True)
                allok &= not diff
    import pickle; pickle.dump(store, open(f'ds3_N{NMAX}.pkl','wb'))
    print('ALL OK' if allok else 'FAILURES')
